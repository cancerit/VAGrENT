package Sanger::CGP::VagrentSV::Base;

use strict;

#----------------------------------------------------------------------------------

 sub open_to_write
 {
     my ($self,$File) = @_;   
     my $WFH;
       
     unless(open($WFH,">$File"))
     {
        print"Failed to open the file:$File because $!\n";
        return(0);
     } 
     
     return($WFH);
}

#----------------------------------------------------------------------------------
 sub open_to_read
 {
     my ($self,$File) = @_;   
     my $rfh;
       
     unless(open($rfh,"$File"))
     {
        print"Failed to open the file:$File because $!\n";
        return(0);
     } 
     
     return($rfh);
}

#----------------------------------------------------------------------------------



sub getSeq {
	my ($self,$fai,$chr,$start,$end)=@_;
	my $dna = $fai->fetch("$chr:$start-$end");
	return $dna;
}



1;