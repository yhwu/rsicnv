#ifndef _LOADDATA_H
#define _LOADDATA_H

void write_rd_to_file(Array<int>& RD, string inpfile);

void concatenate_data(Array<int>& RD);
void patch_N_regions(Array<int>& RD);
void expand_data(Array<int>& RD);
void apply_cap(Array<int>& RD);
void load_data_from_bam(samfile_t *fp_in, bam_index_t *bamidx, int ref,
			Array<int>& RD) ;
void load_data_from_bam(samfile_t *fp_in, bam_index_t *bamidx, 
			int ref, int beg, int end,
			Array<int>& RD) ;
void load_data_from_text(string inpfile, Array<int>& RD) ;

void read_cnvlist(string inpfile, vector<cnv_st>& cnvlist);
void read_cnvlist(string inpfile, vector<cnv_st>& cnvlist, 
		  string & text, vector<size_t>& text_idx);

#endif

