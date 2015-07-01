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
     my ($sel,$File) = @_;   
     my $rfh;
       
     unless(open($rfh,"$File"))
     {
        print"Failed to open the file:$File because $!\n";
        return(0);
     } 
     
     return($rfh);
}

#----------------------------------------------------------------------------------



1;