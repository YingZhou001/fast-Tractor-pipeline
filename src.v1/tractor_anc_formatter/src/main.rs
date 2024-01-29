use std::env;
use std::fmt;
use std::fs;

use flate2::read;
use flate2::write;
use flate2::Compression;
//use std::error::Error;
use std::ffi::OsStr;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Read normal or compressed files seamlessly
/// Uses the presence of a `.gz` extension to decide
pub fn reader(filename: &str) -> Box<dyn BufRead> {
    let path = Path::new(filename);
    let file = match fs::File::open(&path) {
        Err(why) => panic!("couldn't open {}: {:?}", path.display(), why),
        Ok(file) => file,
    };

    if path.extension() == Some(OsStr::new("gz")) {
        Box::new(BufReader::with_capacity(
            128 * 1024,
            read::MultiGzDecoder::new(file),
        ))
    } else {
        Box::new(BufReader::with_capacity(128 * 1024, file))
    }
}

/// Write normal or compressed files seamlessly
/// Uses the presence of a `.gz` extension to decide
// Attempting to have a file writer too
pub fn writer(filename: &str) -> Box<dyn Write> {
    let path = Path::new(filename);
    let file = match fs::File::create(&path) {
        Err(why) => panic!("couldn't open {} : {:?}", path.display(), why),
        Ok(file) => file,
    };

    if path.extension() == Some(OsStr::new("gz")) {
        // Error is here: Created file isn't gzip-compressed
        Box::new(BufWriter::with_capacity(
            128 * 1024,
            write::GzEncoder::new(file, Compression::default()),
        ))
    } else {
        Box::new(BufWriter::with_capacity(128 * 1024, file))
    }
}

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("\tUsage: {} [flare output vcf] [output prefix]\n", &args[0]);
        return Ok(());
        }
    let inpvcf = &args[1];
    let outpref = &args[2];
    let mut anc_num : usize = 0;
    let mut outhap = Vec::new();
    let mut outdosage = Vec::new();
    let mut first_line = Vec::new();

    let reader_file = reader(inpvcf);
    for inpline in reader_file.lines() {
        let mut line = match inpline {
            Ok(line) => line,
                Err(_) => break,
        };
        if &line[0..2] == "##" {
            if &line[0..10] == "##ANCESTRY" {
                anc_num = line[11..].split(',').count();
                println!("ancestry code: {:?}", &line[11..]);
                // creat output file
                for i in 0..anc_num {
                    let mut dfile = String::new();
                    let mut hfile = String::new();
                    fmt::write(&mut dfile, 
                            format_args!("{}.anc{}.dosage.txt.gz", outpref, i))
                        .expect("Error!");
                    fmt::write(&mut hfile, 
                            format_args!("{}.anc{}.hapcount.txt.gz", outpref, i))
                        .expect("Error!");
                    outdosage.push(dfile);
                    outhap.push(hfile);
                }
            }
            continue;
        }
        if &line[0..1] == "#" {
            line = line.replace("#", "");
            line.split('\t').take(5)
                .for_each(|w| first_line.push(w.to_string()));
            line.split('\t').skip(9)
                .for_each(|w| first_line.push(w.to_string()));
            break;
        }
    }
    println!("Dosage files to write: {:?}", outdosage);
    println!("Hapcount files to write: {:?}", outhap);
    
    let mut outhap_filehandle = Vec::new();
    let mut outdosage_filehandle = Vec::new();

    for f in outhap.iter(){
        outhap_filehandle.push(writer(f));
    }

    for f in outdosage.iter(){
        outdosage_filehandle.push(writer(f));
    }

    for i in 0..anc_num {
        let first_line = first_line.join("\t") + "\n";
        outdosage_filehandle[i].write_all(first_line.as_bytes())?;
        outhap_filehandle[i].write_all(first_line.as_bytes())?;
    }

 
    let reader_file = reader(inpvcf);
    for inpline in reader_file.lines() {
        let line = match inpline {
            Ok(line) => line,
            Err(_) => break,
        };
        if &line[0..1] == "#" {continue;}
        let pref_line = line.split('\t').take(5)
            .map(String::from)
            .collect::<Vec<String>>()
            .join("\t") + "\t";
        for i in 0..anc_num {
            let istr = i.to_string();
            outhap_filehandle[i].write_all(pref_line.as_bytes())?;
            outdosage_filehandle[i].write_all(pref_line.as_bytes())?;
            let mut newhap = Vec::new();
            let mut newdosage = Vec::new();
            line.split('\t')
                .skip(9)
                .for_each(|w| {
                    let w_arr = w.split(|ch| ch == '|' || ch == ':')
                    .collect::<Vec<&str>>();
                    let mut h = 0;
                    if w_arr[2] == &istr {h += 1;}
                    if w_arr[3] == &istr {h += 1;}
                    newhap.push(h.to_string());
                    let mut d = 0;
                    if w_arr[2] == &istr && w_arr[0] == "1" {d +=1;}
                    if w_arr[3] == &istr && w_arr[1] == "1" {d +=1;}
                    newdosage.push(d.to_string());
                    });
            outhap_filehandle[i].write_all(newhap.join("\t").as_bytes())?;
            outdosage_filehandle[i].write_all(newdosage.join("\t").as_bytes())?;
            outhap_filehandle[i].write_all(b"\n")?;
            outdosage_filehandle[i].write_all(b"\n")?;
        }
    }        
    
    Ok(())
}
