package Sanger::CGP::VagrentSV::Base;

use strict;

use Log::Log4perl;

use Data::Dumper;
use Capture::Tiny qw(:all);

use FindBin qw($Bin);
Log::Log4perl->init("$Bin/../config/log4perl.vagrentsv.conf");
my $log = Log::Log4perl->get_logger(__PACKAGE__);

use Sanger::CGP::Vagrent qw($VERSION);

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


=head2 _run_cmd
runs external command
Inputs
=over 2
=item cmd - command to run
=back
=cut

sub _run_cmd {
	my($self,$cmd)=@_;
	my ($out,$stderr,$exit)=capture{system($cmd)};
	if($exit) {
			$log->logcroak("Failed to run <<<<<<< \n $cmd  <<<<<< \n with status <<<<<< \n OUT:\n $out  :ERR:\n $stderr EXIT:\n $exit \n <<<<<<< \n");
	}
	else {
		$log->debug("\ncommand <<<<<< \n $cmd \nrun successfully <<<<<<<<< ");
	}
	return $out;
}


1;