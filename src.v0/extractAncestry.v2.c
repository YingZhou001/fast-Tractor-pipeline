/*
 * =====================================================================================
 *
 *       Filename:  extractAncestry.c
 *
 *    Description: Extract ancestry speicifc haplotypes and ancestry dosage from flare.
 *
 *        Version:  1.0
 *        Created:  08/24/2022 09:14:03 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ying Zhou, zy.popg@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/hts.h>

void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    if (strrchr(format, '\n') == NULL) fputc('\n', stderr);
    exit(-1);
}

void print_help(void)
{
  printf("Usage: ./extractAncestry [input: Flare's output vcf.gz] [ancestry code] [output prefix]\n");
  exit(-1);
}


int main(int arg, char **argv)
{
  if(arg == 1 || arg > 5)print_help();

  char *inp=argv[1];
  int code=atoi(argv[2]);
  char *out=argv[3];

  char *outvcf, *outhapanc;
  outvcf = (char*) calloc(strlen(out)+16, sizeof(char));
  outhapanc = (char*) calloc(strlen(out)+16, sizeof(char));

  sprintf(outvcf, "%s.anc%d.vcf.gz", out, code);
  sprintf(outhapanc, "%s.hapanc%d.vcf.gz", out, code);

  htsFile *fp = hts_open(inp, "r");
  if (!fp) error("Failed to open \"%s\" : %s", inp, strerror(errno));
  bcf_hdr_t *hdr = bcf_hdr_read(fp);
  if (!hdr) error("bcf_hdr_read : %s", strerror(errno));
  int sample_num=bcf_hdr_nsamples(hdr);
  int i;
  int r = 0;

  bcf_hdr_t *hdr_out = bcf_hdr_dup(hdr);
  bcf_hdr_remove(hdr_out,BCF_HL_STR,NULL);
  bcf_hdr_remove(hdr_out,BCF_HL_FMT,"AN1");
  bcf_hdr_remove(hdr_out,BCF_HL_FMT,"AN2");

  htsFile *ofp = hts_open(outvcf, "wz");
  if (!ofp) error("Failed to open \"%s\" : %s", outvcf, strerror(errno));
  r = bcf_hdr_write(ofp, hdr_out);
  if ( r != 0 ) error("Failed to write to %s\n", outvcf);

  htsFile *ofp2 = hts_open(outhapanc, "wz");
  if (!ofp2) error("Failed to open \"%s\" : %s", outhapanc, strerror(errno));

  r = bcf_hdr_write(ofp2, hdr_out);
  if ( r != 0 ) error("Failed to write to %s\n", outhapanc);


  int ann1 = 0; 
  int ann2 = 0; 
  int gtn = 0; 
  int tmp = 0;
  int32_t *an1 = NULL;
  int32_t *an2 = NULL;
  int32_t *gt = NULL;
  int32_t *tmpia = (int*)malloc(sample_num*2*sizeof(int));

  bcf1_t *lbuf =  bcf_init();

  while(1 > 0){
    r = bcf_read(fp, hdr, lbuf);
    if(r != 0) break;
    bcf_get_format_int32(hdr, lbuf, "AN1", &an1, &ann1);
    bcf_get_format_int32(hdr, lbuf, "AN2", &an2, &ann2);
    bcf_get_genotypes(hdr, lbuf, &gt, &gtn);

    for (i=0; i< sample_num; i++) {
      if(an1[i] == code)tmpia[2*i]=gt[2*i]; 
      else {
	tmpia[2*i]=bcf_gt_missing|1;    
      }
      if(an2[i] == code)tmpia[2*i+1]=gt[2*i+1];
      else {
	tmpia[2*i+1]=bcf_gt_missing|1;    
      }
    }
    bcf_update_genotypes(hdr_out, lbuf, tmpia, sample_num*2);
    bcf_update_format_int32(hdr_out, lbuf, "AN1", NULL, 0);
    bcf_update_format_int32(hdr_out, lbuf, "AN2", NULL, 0);

    r = bcf_write1(ofp, hdr_out, lbuf);
    //    printf("%ld %d\n", lbuf->pos, r);

    for (i=0; i< sample_num; i++) {
      if(an1[i] == code)tmpia[2*i]=bcf_gt_phased(1); 
      else tmpia[2*i]=bcf_gt_phased(0);
      if(an2[i] == code)tmpia[2*i+1]=bcf_gt_phased(1); 
      else tmpia[2*i+1]=bcf_gt_phased(0);
    }

    bcf_update_genotypes(hdr_out, lbuf, tmpia, sample_num*2);
    r = bcf_write1(ofp2, hdr_out, lbuf);
    bcf_clear(lbuf);

  }

  hts_close(ofp);
  hts_close(ofp2);
  bcf_destroy(lbuf);

  free(tmpia);
  free(an1);
  free(an2);
  free(gt);
  free(outvcf);
  free(outhapanc);

  bcf_hdr_destroy(hdr);
  bcf_hdr_destroy(hdr_out);
  hts_close(fp);

  return 0;
}
