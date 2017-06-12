#!/usr/bin/perl
BEGIN {
  $SIG{__WARN__} = sub {warn $_[0] unless(( $_[0] =~ m/^Subroutine Tabix.* redefined/) || ($_[0] =~ m/^Use of uninitialized value \$buf/))};
};

use strict;
use Const::Fast qw(const);
use Tabix;
use Data::Dumper;


const my $FLANK_DIST => 750000;#750kb
const my $BED_FORM => qq{%s\t%d\t%d\t%s\t.\t%s\t%s\t%s\t%s\n}; # chr, start, stop, name, (score), strand, other, chrArmi, centromerLoc


my ($gencode_gff3, $enhancer_bed, $fai, $cytoband, $rem_chr_prefix, $chr_restrict) = @ARGV;

generate_data_files($gencode_gff3, $enhancer_bed, $fai, $cytoband, $rem_chr_prefix, $chr_restrict);

sub generate_data_files {
  my ($gencode_gff3, $enhancer_bed,$fai,$cytoband,$rem_chr_prefix, $chr_restrict) = @_;
  my $chr_len = _getHash($fai);
  my $cyto_len = _getHash($cytoband);
  print Dumper $cyto_len;
  my $gene_enhancer_bed='gene_enhancer.bed';
	open my $gene_enhancer_fh, '>', $gene_enhancer_bed;
	_parse_enhancer($enhancer_bed, $rem_chr_prefix, $chr_restrict, $chr_len,$gene_enhancer_fh,$cyto_len);
  _parse_gene($gencode_gff3, $rem_chr_prefix, $chr_restrict, $chr_len,$gene_enhancer_fh,$cyto_len);
  # sorts input file
	_sort_file([$gene_enhancer_bed]);

}

sub _parse_gene {
  my ($gff3, $rem_chr_prefix, $chr_restrict, $chr_len,$gene_enhancer_fh,$cyto_len) = @_;
  open my $IN, '<', $gff3 or die "Failed to open $gff3: $!";
  my $gene_id = 1;
  while(<$IN>) {
    chomp $_;
    my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $attr) = split /\t/, $_;
    $chr =~ s/^chr// if($rem_chr_prefix == 1);
    next if(defined $chr_restrict && $chr ne $chr_restrict);
		my $startL=($start-$FLANK_DIST)-1;
		my $endH=$end+$FLANK_DIST;
    # upstream flank of $FLANK_DIST
    if($start < $FLANK_DIST){$startL=0;}
		if($chr_len->{$chr} < $endH) {
			$endH=$chr_len->{$chr};
		}
		my $centromere_loc=$cyto_len->{$chr};
		# upstream flank of $FLANK_DIST
		my($chr_arm)=_getChrArm($cyto_len->{$chr},$start);
    printf $gene_enhancer_fh $BED_FORM, $chr, $startL, $start-1, '-g'.$gene_id, $strand, $attr, $chr_arm, $centromere_loc;
    # the original gene coords
    printf $gene_enhancer_fh $BED_FORM, $chr, $start-1, $end, 'g'.$gene_id, $strand, $attr, $chr_arm, $centromere_loc;
    # downstream flank of $FLANK_DIST
    printf $gene_enhancer_fh $BED_FORM, $chr, $end, $endH, '+g'.$gene_id, $strand, $attr, $chr_arm, $centromere_loc;
    $gene_id++;
  }
  close $IN;
}

sub _parse_enhancer {
  my ($bed, $rem_chr_prefix, $chr_restrict, $chr_len,$gene_enhancer_fh,$cyto_len) = @_;
  open my $IN, '<', $bed or die "Failed to open $bed: $!";
  my $enhancer_id = 1;
  while(<$IN>) {
    chomp $_;
    my ($chr, $start, $end, $relation) = split /\t/, $_;
    $chr =~ s/^chr// if($rem_chr_prefix == 1);

    next if(defined $chr_restrict && $chr ne $chr_restrict);
    # the core enhancer location is irrelevant
    my $startL=($start-$FLANK_DIST)-1;
		my $endH=$end+$FLANK_DIST;
    # upstream flank of $FLANK_DIST
    if($start < $FLANK_DIST){$startL=0;}
		if($chr_len->{$chr} < $endH) {
			$endH=$chr_len->{$chr};
		}
		my $centromere_loc=$cyto_len->{$chr};
		# upstream flank of $FLANK_DIST
		my($chr_arm)=_getChrArm($cyto_len->{$chr},$start);
    printf $gene_enhancer_fh $BED_FORM, $chr, $startL, $start, '-e'.$enhancer_id, '.', $relation, $chr_arm, $centromere_loc;
    # the original enhancer coords
    printf $gene_enhancer_fh $BED_FORM, $chr, $start-1, $end, 'e'.$enhancer_id, '.', $relation, $chr_arm,$centromere_loc;
    # downstream flank of $FLANK_DIST
    printf $gene_enhancer_fh $BED_FORM, $chr, $end, $endH, '+e'.$enhancer_id, '.', $relation, $chr_arm,$centromere_loc;
    $enhancer_id++;
  }
  close $IN;
}

sub _sort_file {
	my($files)=@_;
	foreach my $file (@$files) {
		system("sort -k1,1 -k2,2n $file  >sorted_$file");
		}
}

# ---------------------
sub _getHash {
  my($file_name)=@_;
  my $chr_length;
  open my $fh, '<', $file_name or die 'unable to open file'.$!;
  while(<$fh>){
    my($key,$val)=split "\t";
    $key=~s/chr//;
    $chr_length->{$key}=$val;
  }
 $chr_length;
}


sub _getChrArm{
 my($cytoband_pos,$genePos)=@_;
	my $gene_arm='p';
	if ($cytoband_pos < $genePos){
		$gene_arm='q';
	}
	$gene_arm;
}

__DATA__

perl ~/sv_annot.pl /var/tmp/gencode19_protein_coding_gene.gff3 /var/tmp/roadmap_stringent_enhancers_with_Ensembl_ID_and_name.bed genome_h37.fa.fai centrometerPos.txt 1 21
