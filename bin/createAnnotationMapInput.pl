#!/usr/bin/perl
BEGIN {
  $SIG{__WARN__} = sub {warn $_[0] unless(( $_[0] =~ m/^Subroutine Tabix.* redefined/) || ($_[0] =~ m/^Use of uninitialized value \$buf/))};
};

use strict;
use Const::Fast qw(const);
use Tabix;

const my $FLANK_DIST => 750000;#750kb
const my $BED_FORM => qq{%s\t%d\t%d\t%s\t.\t%s\t%s\n}; # chr, start, stop, name, (score), strand, other


my ($gencode_gff3, $enhancer_bed, $fai, $rem_chr_prefix, $chr_restrict) = @ARGV;

generate_data_files($gencode_gff3, $enhancer_bed, $fai, $rem_chr_prefix, $chr_restrict);

sub generate_data_files {
  my ($gencode_gff3, $enhancer_bed, $fai, $rem_chr_prefix, $chr_restrict) = @_;
  my $chr_len = _getHash($fai);
  my $gene_enhancer_bed='gene_enhancer.bed';
	open my $gene_enhancer_fh, '>', $gene_enhancer_bed;
	_parse_enhancer($enhancer_bed, $rem_chr_prefix, $chr_restrict, $chr_len,$gene_enhancer_fh);
  _parse_gene($gencode_gff3, $rem_chr_prefix, $chr_restrict, $chr_len,$gene_enhancer_fh);
  # sorts input file
	_sort_file([$gene_enhancer_bed]);

}

sub _parse_gene {
  my ($gff3, $rem_chr_prefix, $chr_restrict, $chr_len,$gene_enhancer_fh) = @_;
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
		# upstream flank of $FLANK_DIST
    printf $gene_enhancer_fh $BED_FORM, $chr, $startL, $start-1, '-g'.$gene_id, $strand, $attr;
    # the original gene coords
    printf $gene_enhancer_fh $BED_FORM, $chr, $start-1, $end, 'g'.$gene_id, $strand, $attr;
    # downstream flank of $FLANK_DIST
    printf $gene_enhancer_fh $BED_FORM, $chr, $end, $endH, '+g'.$gene_id, $strand, $attr;
    $gene_id++;
  }
  close $IN;
}

sub _parse_enhancer {
  my ($bed, $rem_chr_prefix, $chr_restrict, $chr_len,$gene_enhancer_fh) = @_;
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
		# upstream flank of $FLANK_DIST
    printf $gene_enhancer_fh $BED_FORM, $chr, $startL, $start, '-e'.$enhancer_id, '.', $relation;
    # the original enhancer coords
    printf $gene_enhancer_fh $BED_FORM, $chr, $start-1, $end, 'e'.$enhancer_id, '.', $relation;
    # downstream flank of $FLANK_DIST
    printf $gene_enhancer_fh $BED_FORM, $chr, $end, $endH, '+e'.$enhancer_id, '.', $relation;
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
    $chr_length->{$key}=$val;
  }
 $chr_length;
}



__DATA__

perl ~/sv_annot.pl /var/tmp/gencode19_protein_coding_gene.gff3 /var/tmp/roadmap_stringent_enhancers_with_Ensembl_ID_and_name.bed genome_h37.fa.fai 1 21
