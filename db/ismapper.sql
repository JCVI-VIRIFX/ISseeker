
##
## Old unused tables
##
DROP TABLE IF EXISTS SNP;
DROP TABLE IF EXISTS corepos;
DROP TABLE IF EXISTS insertion_seqs;
DROP TABLE IF EXISTS facX;
DROP TABLE IF EXISTS plasmid_matches;
DROP TABLE IF EXISTS contigs;
DROP TABLE IF EXISTS match_aggregates;
DROP TABLE IF EXISTS plasmid_matches;


DROP TABLE IF EXISTS is_flank;
DROP TABLE IF EXISTS is_query_feature;
DROP TABLE IF EXISTS is_run;

CREATE TABLE is_run
(
        is_run_id int NOT NULL AUTO_INCREMENT,
	is_run_element VARCHAR(255) NOT NULL,
	is_run_genome VARCHAR(255) NOT NULL,
	is_run_references VARCHAR(512) NOT NULL,
	is_req_pct_id float NOT NULL,
	flank_req_pct_id float NOT NULL,
	PRIMARY KEY(is_run_id),
	UNIQUE INDEX (is_run_element, is_run_genome, is_run_references)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


##
## is_query_feature - General BLAST information
## 
## Types of entries in this table:
##
##   flank blasted against annotation (type=F)
##   IS blasted against annotation (type=E)
##  
##  Not saved:
##    IS blasted against contigs to generate flanks
 
CREATE TABLE is_query_feature
(
        is_qf_id int NOT NULL AUTO_INCREMENT,
        is_run_id int NOT NULL,
        is_element VARCHAR(255) NOT NULL,
        reference VARCHAR(255) NOT NULL,
        q_genome VARCHAR(255) NOT NULL,
        s_genome VARCHAR(255) NOT NULL,
        is_pct_id float NOT NULL,      -- pct ID of IS vs. contig
        flank_pct_id float NULL,       -- pct ID of the flank vs. annotation
        contig_name VARCHAR(255) NOT NULL,
        feat_type ENUM('F', 'E','X')  NOT NULL comment 'Flank/Embedded IS/Unknown',
        flank_begin_end ENUM('B','E','X')  NOT NULL comment  'Begin/End/Unknown',
        orientation ENUM('F','R','X') NOT NULL comment 'Forward/Reverse/Unknown',
        begin_base int NOT NULL,
        end_base int NOT NULL,
        is_annotated bit NOT NULL comment 'TRUE = has a record in the is_flank table',
        PRIMARY KEY(is_qf_id),
        INDEX(begin_base,end_base),
        UNIQUE INDEX (is_run_id, is_element, reference, q_genome, s_genome, contig_name, begin_base )
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

ALTER TABLE is_query_feature ADD CONSTRAINT qf_run_fk FOREIGN KEY (is_run_id) REFERENCES is_run(is_run_id);



## 
##  contig_name moved to is_query_feature
##
CREATE TABLE is_flank
(
    is_flank_id int NOT NULL AUTO_INCREMENT,
    is_run_id int NOT NULL,
    is_qf_id int NOT NULL,
    is_element VARCHAR(255) NOT NULL,
    reference VARCHAR(255) NOT NULL,
    q_genome VARCHAR(255) NOT NULL,
    flank_pct_id float NOT NULL,
    begin_end ENUM('B','E','X')  NOT NULL comment 'Begin/End/Unknown',
    orientation ENUM('F','R','X') NOT NULL comment 'Forward/Reverse/Unknown',
    nearest_base int NOT NULL,
    after_gene VARCHAR(255) NULL,
    in_gene VARCHAR(255) NULL,
    before_gene VARCHAR(255) NULL,
    match_quality ENUM('W', 'P', 'X') NOT NULL comment 'Whole/Partial/None',
    flank_separation int NULL,
    mate_flank_id int NULL,
    contig_count int NOT NULL comment 'Number of contigs that matched this flank (only 1 is in is_query_feature)',
    PRIMARY KEY(is_flank_id),
    INDEX(reference,nearest_base),
    UNIQUE INDEX (is_element, reference, q_genome, nearest_base )
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

ALTER TABLE is_flank ADD CONSTRAINT fk_run_fk FOREIGN KEY (is_run_id) REFERENCES is_run(is_run_id);
ALTER TABLE is_flank ADD CONSTRAINT fk_mate_flank_fk FOREIGN KEY (mate_flank_id) REFERENCES is_flank(is_flank_id);
ALTER TABLE is_flank ADD CONSTRAINT fk_query_feat_fk FOREIGN KEY (is_qf_id) REFERENCES is_query_feature(is_qf_id);



CREATE OR REPLACE VIEW flank_pairs AS

SELECT a.orientation as orientation, 
        a.is_flank_id as begin_flank_id, 
        a.is_element as is_element, 
        a.q_genome as q_genome,
        a.reference as reference, 
        a.nearest_base as begin_base,  
        a.after_gene as begin_after_gene, 
        a.in_gene as begin_in_gene, 
        a.before_gene as begin_before_gene,       
        a.flank_separation, 
        a.flank_pct_id as begin_flank_pct_id,
        a.match_quality as begin_match_quality,       
        b.nearest_base as end_base, 
        b.after_gene as end_after_gene, 
        b.in_gene as end_in_gene, 
        b.before_gene as end_before_gene, 
        b.flank_pct_id as end_flank_pct_id,
        b.match_quality as end_match_quality, 
        b.is_flank_id as end_flank_id,
		least(a.nearest_base,b.nearest_base) AS lowest_base
from is_flank a JOIN is_flank b on a.mate_flank_id = b.is_flank_id
WHERE a.begin_end = 'B'
UNION
SELECT a.orientation as orientation,
        a.is_flank_id as begin_flank_id, 
        a.is_element as is_element, 
        a.q_genome as q_genome,
        a.reference as reference, 
        a.nearest_base as begin_base, 
        a.after_gene as begin_after_gene, 
        a.in_gene as begin_in_gene, 
        a.before_gene as begin_before_gene,        
        a.flank_separation, 
        a.flank_pct_id as begin_flank_pct_id,
        a.match_quality as begin_match_quality, 
        b.nearest_base as end_base, 
        b.after_gene as end_after_gene,
        b.in_gene as end_in_gene, 
        b.before_gene as end_before_gene, 
        b.flank_pct_id as end_flank_pct_id,
        b.match_quality as end_match_quality, 
        b.is_flank_id as end_flank_id,
		least(COALESCE(a.nearest_base,b.nearest_base),COALESCE(b.nearest_base,a.nearest_base)) AS lowest_base
from is_flank a  LEFT OUTER JOIN is_flank b on a.mate_flank_id = b.is_flank_id where a.begin_end = 'B' and a.mate_flank_id is null
UNION
select b.orientation as orientation,
        a.is_flank_id as begin_flank_id, 
        b.is_element as is_element, 
        b.q_genome as q_genome,
        b.reference as reference, 
        a.nearest_base as begin_base, 
        a.after_gene as begin_after_gene, 
        a.in_gene as begin_in_gene, 
        a.before_gene as begin_before_gene,         
        a.flank_separation, 
        a.flank_pct_id as begin_flank_pct_id,
        a.match_quality as begin_match_quality, 
        b.nearest_base as end_base, 
        b.after_gene as end_after_gene,
        b.in_gene as end_in_gene, 
        b.before_gene as end_before_gene, 
        b.flank_pct_id as end_flank_pct_id,
        b.match_quality as end_match_quality, 
        b.is_flank_id as end_flank_id,
		least(COALESCE(a.nearest_base,b.nearest_base),COALESCE(b.nearest_base,a.nearest_base)) AS lowest_base
from is_flank a  RIGHT OUTER JOIN is_flank b on a.mate_flank_id = b.is_flank_id where b.begin_end = 'E' and b.mate_flank_id is null ;
