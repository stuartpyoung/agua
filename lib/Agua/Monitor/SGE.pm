use MooseX::Declare;

#### DEBUG

=head2

	PACKAGE		Agua::Monitor::SGE

    VERSION:        0.03

    PURPOSE

        1. MONITOR JOBS RUN ON AN SGE (PORTABLE BATCH SCHEDULER) SYSTEM

	HISTORY

		0.03	Streamlined by removing SQLite-specific methods
		0.02	Added monitor_jobs failback on 'ERROR' qstat output
		0.01	Basic version

=cut 

use strict;
use warnings;
use Carp;

use FindBin qw($Bin);
use lib "$Bin/../..";

#### EXTERNAL MODULES
use POSIX;
use Data::Dumper;
use DBI;
use File::Path;

#### INTERNAL MODULES
#class Agua::Monitor::SGE with Agua::Common::SGE {
class Agua::Monitor::SGE with Agua::Common::SGE {

#### EXTERNAL MODULES
use Data::Dumper;

# INTS
has 'pid'		=> ( isa => 'Int|Undef', is => 'rw' );
has 'tries'		=> ( isa => 'Int|Undef', is => 'rw', default => 20 ); 	#### BEFORE QUIT
has 'sleep'		=> ( isa => 'Int|Undef', is => 'rw', default => 5 );
# STRINGS
has 'cluster'	=> ( isa => 'Str|Undef', is => 'rw' );
has 'username'  => ( isa => 'Str|Undef', is => 'rw' );
has 'workflow'  => ( isa => 'Str|Undef', is => 'rw' );
has 'project'   => ( isa => 'Str|Undef', is => 'rw' );
has 'qsub'		=> ( isa => 'Str|Undef', is => 'rw', required => 0 );
has 'qstat'		=> ( isa => 'Str|Undef', is => 'rw', required => 0 );
has 'errorregex'=> ( isa => 'Str|Undef', is => 'rw', default => "(error|unable to contact)");
has 'jobidregex'=> ( isa => 'Str|Undef', is => 'rw', default => qq{^(Your job|Your job-array) (\\S+)});
# OBJECTS
has 'conf'		=> ( isa => 'Conf::Agua', is => 'rw', required => 1 );
has 'dbobject'	=> ( isa => 'Agua::DBase::MySQL', is => 'rw', required => 0 );

#####///}}

method BUILD ($hash) {

	##### SET DATABASE HANDLE
	#$self->setEnv();

	#### SET qstat EXECUTABLE LOCATION
	$self->setQstat();
}


method setQstat {
	my $conf = $self->conf();

	my $qstat = $conf->getKeyValue("cluster", "QSTAT");
	$qstat = $conf->getKeyValue("cluster", "QSTAT") if not defined $qstat;
	print "Agua::Monitor::SGE::setEnv    sgeroot not defined\n" and exit if not defined $qstat;

	$self->qstat($qstat);	
}

=head2

	SUBROUTINE		setEnv

	PURPOSE

		SET SGE ENVIRONMENT VARIABLES

=cut

method setEnv {


	my $conf = $self->conf();

	my $sgeroot = $conf->getKeyValue('agua', "SGEROOT");
	$sgeroot = $conf->getKeyValue("cluster", "SGEROOT") if not defined $sgeroot;
	my $qmasterport = $conf->getKeyValue("cluster", "SGEQMASTERPORT");
	$qmasterport = $conf->getKeyValue("cluster", "SGEQMASTERPORT") if not defined $qmasterport;
	my $execdport = $conf->getKeyValue("cluster", "SGEEXECDPORT");
	$execdport = $conf->getKeyValue("cluster", "SGEEXECDPORT") if not defined $execdport;

	#### CHECK INPUTS
	print "Agua::Monitor::SGE::setEnv    sgeroot not defined\n" and exit if not defined $sgeroot;
	print "Agua::Monitor::SGE::setEnv    qmasterport not defined\n" and exit if not defined $qmasterport;
	print "Agua::Monitor::SGE::setEnv    execdport not defined\n" and exit if not defined $execdport;

	$ENV{'SGE_ROOT'} = $sgeroot;
	$ENV{'SGE_QMASTER_PORT'} = $qmasterport;
	$ENV{'SGE_EXECD_PORT'} = $execdport;

}

=head2

	SUBROUTINE		submitJob

	PURPOSE

		SUBMIT A JOB TO THE CLUSTER AND RETURN THE JOB ID

		IF JOB ID IS UNDEFINED, REPORT AND QUIT

=cut
method submitJob ($args) {

#exit;


	my $queue 		= 	$args->{queue};
	my $batch 		= 	$args->{batch};
	$batch = '' if not defined $batch;
	my $scriptfile 	= 	$args->{scriptfile};
	my $qsub 		=	$args->{qsub};

	#### CHECK INPUTS
	print "Agua::Monitor::SGE::submitJob   queue not defined. Returning\n" and return if not defined $queue;
	print "Agua::Monitor::SGE::submitJob   qsub not defined. Returning\n" and return if not defined $qsub;
	print "Agua::Monitor::SGE::submitJob   Scriptfile not defined. Returning null\n" and return if not defined $args->{scriptfile};

	#### SET SGE_CELL TO DETERMINE WHICH CLUSTER TO RUN JOB ON
	my $cluster = $args->{cluster};
	$ENV{'SGE_CELL'} = $cluster if defined $cluster;

	#### SET COMMAND
	my $command;

	#### SET ENVIRONMENT VARIABLES
	$command .= $args->{envars}->{tostring} . ";" if defined $args->{envars};
	if ( not defined $args->{envars} )
	{
		$command .= "export SGE_CELL=" . $ENV{'SGE_CELL'} . ";" if defined $ENV{'SGE_CELL'}; 
		$command .= "export SGE_ROOT=" . $ENV{'SGE_ROOT'} . ";" if defined $ENV{'SGE_ROOT'}; 	$command .= "export SGE_QMASTER_PORT=" . $ENV{'SGE_QMASTER_PORT'} . ";" if defined $ENV{'SGE_QMASTER_PORT'}; 
		$command .= "export SGE_EXECD_PORT=" . $ENV{'SGE_EXECD_PORT'} . ";" if defined $ENV{'SGE_EXECD_PORT'}; 
	}

	#### SET QSUB LINE
	$command .= "$qsub $batch -q $queue $scriptfile";



#return "7.1-1:1";



	#### SUBMIT
	my $output = `$command`;


	##### LATER: POSSIBLY REPLACE WITH OPEN PIPE METHOD (FROM CPAN MODULE Schedule::SGE)
	#open(QSUB, "|$pipe > /tmp/$$.out 2>&1") || die "Can't open the pipe to submit jobs to";
	#close QSUB;
	#return 0 unless (-e "/tmp/$$.out"); 
	#open (IN, "/tmp/$$.out") || die "Can't open /tmp/$$.out";
	#my $line=<IN>;
	#close IN;

	#### SUBMIT COMMAND OUTPUT FORMAT:
	#### Your job-array 214.1-10:1 ("bowtie-chr22") has been submitted
	my $jobid_regex = $self->jobidregex();

	use re 'eval';	# EVALUATE AS REGEX
	my ($job_type, $job_id) = $output =~ /$jobid_regex/;
	no re 'eval';	# STOP EVALUATE AS REGEX
	print " { error: 'Agua::Monitor::SGE::submitJob    job_id not defined' }"
		and exit if not defined $job_id;


	return $job_id;	
}


=head2

    SUBROUTINE      remainingJobs

    PURPOSE

        RETURN THE LIST OF PIDS CURRENTLY IN THE QUEUE

	NOTES

		FINISHED JOBS DISAPPEAR IMMEDIATELY FROM THE qstat

		OUTPUT BUT THERE MAY BE A LAG BETWEEN SUBMISSION AND

		THE JOB APPEARING IN THE qstat OUTPUT SO TRY SEVERAL

		TIMES IF A JOB DOESN'T APPEAR TO BE IN THE QUEUE.

=cut

method remainingJobs ($job_ids) {


	return $self->remainingBatchJobs($job_ids) if $$job_ids[0] =~ /^\d+\.\d+/;


	##### GET LIST OF JOBS IN QSTAT
	my $matched = [];
	my $tries = 3;
	my $lines = $self->jobLines();
	while ( @$job_ids and $tries )
	{
		for ( my $i = 0; $i < @$job_ids; $i++ )
		{
			my $job_id = $$job_ids[$i];
			foreach my $line ( @$lines )
			{
				print "Agua::Monitor::SGE::remainingJobs    line: $line\n";
				if ( $line =~ /^\s*$job_id\s+/ )
				{
					print "Agua::Monitor::SGE::remainingJobs    matched\n";
					push @$matched, splice @$job_ids, $i, 1;
					$i--;
					last;
				}
			}

			print "Agua::Monitor::SGE::remainingJobs    job_id $i matched: @$matched\n";
		}

		last if not @$job_ids;

		$tries--;
		sleep(5);
	}

    return $matched;
}


=head2

    SUBROUTINE      remainingBatchJobs

    PURPOSE

        RETURN THE LIST OF PIDS CURRENTLY IN THE QUEUE

	NOTES

		PARSE OUT JOBARRAY IDS FROM qstat OUTPUT TO DETERMINE

		IF JOBS ARE STILL QUEUED OR RUNNING. RETURN THE LIST 

		OF JOBS STILL QUEUED OR RUNNING

=cut

method remainingBatchJobs ($job_ids){	

	my $matched = [];
	my $tries = 3;
	my $lines = $self->jobLines();


	while ( @$job_ids and $tries )
	{

		for ( my $i = 0; $i < @$job_ids; $i++ )
		{
			my ($job_id, $task_info) = $$job_ids[$i] =~ /^(\d+)\.(.+)$/;

			my $tasks = $self->getTasks($task_info);

			#### JOB ID FORMAT:
			#### Your job-array 214.1-10:1 ("bowtie-chr22") has been submitted
			my $hit = 0;
			foreach my $line ( @$lines )
			{
				next if not $line =~ /^\s*(\d+)\s+.+?(\S+)\s*$/;
				my $current_job_id = $1;
				my $current_task_info = $2;

				if ( $job_id == $current_job_id )
				{

					my $current_tasks = $self->getTasks($current_task_info);

					foreach my $task ( @$tasks )
					{
						foreach my $current_task (@$current_tasks)
						{
							if ( $task == $current_task )
							{
								$hit = 1;
								last;
							}
						}

						last if $hit;

					}	#### COMPARE TASKS

					if ( $hit )
					{

						push @$matched, splice @$job_ids, $i, 1;
						$i--;
						last;
					}

				}	#### COMPARE job_id AND current_job_id

				last if $hit;

			}	#### lines			

		}

		last if not @$job_ids;

		$tries--;
		sleep(5);
	}


    return $matched;
}


method getTasks ($task_info) {

	my $tasks;
	push @$tasks, $task_info and return $tasks if $task_info =~ /^\d+$/;

	push @$tasks, split ",", $task_info and return $tasks if $task_info =~ /^\d+,\d+$/;

	#### HANDLE BATCH JOB
	my ($task_ids, $step) = $task_info =~ /^(\S+)\:(\d+)$/;
	if ( $task_ids =~ /^(\d+)-(\d+)$/ )
	{
		my $start = $1;
		my $stop = $2;

		for ( my $i = $start; $i <= $stop; $i+=$step )
		{
			push @$tasks, $i;
		}
	}


	return $tasks;
}



method jobIds {


	#### GET LIST OF JOBS IN QSTAT
	my $lines = $self->jobLines();

	#### PARSE OUT IDS FROM LIST
	my $job_ids = [];
	foreach my $line ( @$lines )
	{
		use re 'eval';# EVALUATE AS REGEX
		my $jobid_regex = $self->jobidregex();
		$line =~ /$jobid_regex/;
		push @$job_ids, $1 if defined $1 and $1;
		no re 'eval';# STOP EVALUATING AS REGEX
	}

	return $job_ids;
}


=head2

	SUBROUTINE		jobLines

	PURPOSE

		RETURN THE LINES FROM A QSTAT CALL:

			EXEC QSTAT COMMAND AND COLLECT RESULT LINES AS THEY <STREAM> OUT

=cut

method jobLines {


	#### TRY REPEATEDLY TO GET A CLEAN QSTAT REPORT
	my $sleep		=	$self->sleep();	
	my $tries		=	$self->tries();	

	#### 
	#### GET ERROR MESSAGES ALSO
	my $sgebin = $self->sgebinCommand();
	my $command = "$sgebin/qstat 2>&1 |";

	#### REPEATEDLY TRY SYSTEM CALL UNTIL A NON-ERROR RESPONSE IS RECEIVED
	my $error_regex = $self->errorregex();

	#my $result = $self->repeatTries($command, $sleep, $tries);

	my $result = qq{job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID 
-----------------------------------------------------------------------------------------------------------------
     95 0.00000 test       root         qw    05/17/2011 05:33:18                                    1        
     96 0.00000 test       root         qw    05/17/2011 05:33:18                                    1        
     97 0.00000 test       root         qw    05/17/2011 05:33:18                                    1        
     98 0.00000 test       root         qw    05/17/2011 05:33:18                                    1        
     99 0.00000 test       root         qw    05/17/2011 05:33:18                                    1        
    100 0.00000 test       root         qw    05/17/2011 05:33:18                                    1        
    101 0.00000 test       root         qw    05/17/2011 05:33:18                                    1        
    102 0.00000 test       root         qw    05/17/2011 05:33:18                                    1        
    103 0.00000 test       root         qw    05/17/2011 05:33:18                                    1        
    104 0.00000 test       root         qw    05/17/2011 05:33:18                                    1        
    105 0.00000 test       root         qw    05/17/2011 05:33:18                                    1        
    106 0.00000 test       root         qw    05/17/2011 05:33:18                                    1        
    107 0.00000 test       root         qw    05/17/2011 05:33:18                                    1        
    108 0.00000 test       root         qw    05/17/2011 05:33:18                                    1        
    109 0.00000 test       root         qw    05/17/2011 05:33:19                                    1        
    110 0.00000 test       root         qw    05/17/2011 05:33:19                                    1        
    111 0.00000 test       root         qw    05/17/2011 05:33:19                                    1        
    112 0.00000 test       root         qw    05/17/2011 05:33:19                                    1        
    113 0.00000 test       root         qw    05/17/2011 05:33:19                                    1        
    114 0.00000 test       root         qw    05/17/2011 05:33:19                                    1        
    115 0.00000 test       root         qw    05/17/2011 05:33:19                                    1        
    116 0.00000 test       root         qw    05/17/2011 05:33:19                                    1        
    117 0.00000 test       root         qw    05/17/2011 05:33:19                                    1        
    118 0.00000 test       root         qw    05/17/2011 05:33:19                                    1        
    119 0.00000 test       root         qw    05/17/2011 05:33:19                                    1        
    120 0.00000 test       root         qw    05/17/2011 05:33:19                                    1        
    121 0.00000 test       root         qw    05/17/2011 05:33:19                                    1        
    122 0.00000 test       root         qw    05/17/2011 05:33:19                                    1        
    123 0.00000 test       root         qw    05/17/2011 05:33:19                                    1
};

	my @lines = split "\n", $result;
	my $jobs_list = [];
	foreach my $line ( @lines )
	{
		use re 'eval';# EVALUATE AS REGEX
		my $jobid_regex = $self->jobidregex();
		push @$jobs_list, $line if $line =~ /$jobid_regex/;
		no re 'eval';# EVALUATE AS REGEX
	}

	return \@lines;
}


=head2

	SUBROUTINE		repeatTries

	PURPOSE

		REPEATEDLY TRY SYSTEM CALL UNTIL A NON-ERROR RESPONSE IS RECEIVED

=cut

method repeatTries ($command, $sleep, $tries) {


	my $result = '';	
	my $error_message = 1;
	while ( $error_message and $tries )
	{
		open(COMMAND, $command) or die "Can't exec command: $command\n";
		while(<COMMAND>) {
			$result .= $_;
		}
		close (COMMAND);

		use re 'eval';	# EVALUATE AS REGEX
		my $error_regex = $self->errorregex();
		$error_message = $result =~ /$error_regex/;
		no re 'eval';# STOP EVALUATING AS REGEX

		#### DECREMENT TRIES AND SLEEP
		$tries--;
		print ".";


		sleep($sleep) if $error_message;
	}

	return $result;
}


=head2

    SUBROUTINE      jobLineStatus

    PURPOSE

        RETURN THE STATUS OF A PARTICULAR JOB IDENTIFIED BY JOB ID

	NOTES

		CALLED IN Sampler.pm BY fastaInfos AND printFiles SUBROUTINES

=cut

method jobLineStatus {
	my $job_id    	=	shift;
    my $job_lines	=   shift;



	#### GET LIST OF JOBS IN QSTAT
	my $lines = $self->jobLines();

	my $status;
    foreach my $line ( @$lines )
    {
		my $job_id_match = $line =~ /^$job_id\./;
		if ( $job_id_match )
		{
			$status = $self->lineStatus($line);
			last;
        }
    }

    return $status;
}



=head2

    SUBROUTINE      jobStatus

    PURPOSE

        RETURN THE STATUS OF A PARTICULAR JOB IDENTIFIED BY JOB ID

	NOTES

		CALLED IN Sampler.pm BY fastaInfos AND printFiles SUBROUTINES

=cut

method jobStatus ($job_id) {

	my $jobcheck = $self->checkJob($job_id);	

	my $status = $self->checkToStatus($jobcheck);

    return $status;
}

=head2

	SUBROUTINE		checkJob

	PURPOSE

		GET RESPONSE FROM checkJob CALL


=cut

method checkJob ($job_id) {

	die "Agua::Monitor::SGE::job_id    job_id not defined. Exiting.\n" if not defined $job_id;

	#qstat -a, -i, -r, -u, -n, -s, -G or -M 
	my $qstat = $self->qstat();
	die "Agua::Monitor::SGE::qstat    qstat not defined. Exiting.\n" if not defined $qstat;

	#### SET COMMAND
	my $command = "$qstat -j $job_id 2>&1 |";

	#### SLEEP BETWEEN TRIES
	my $sleep		=	$self->sleep();	
	my $tries		=	$self->tries();	

	#### REPEATEDLY TRY SYSTEM CALL UNTIL A NON-ERROR RESPONSE IS RECEIVED
	#### OR UNTIL A PRESET NUMBER OF TRIES HAVE ELAPSED
	my $result = $self->repeatTries($command, $sleep, $tries);

	return $result;	
}

=head2

	SUBROUTINE		checkToStatus

	PURPOSE

		PARSE qstat OUTPUT TO GET JOB STATUS


=cut

method checkToStatus {
	my $jobcheck	=	shift;



	my $statusHash;
	($statusHash->{status}) 	= $jobcheck =~ /State:\s*(\S+)/ms;
	($statusHash->{submittime}) = $jobcheck =~ /submission_time:\s*([^\n]+)/ms;
	($statusHash->{starttime}) 	= $jobcheck =~ /StartTime:\s*([^\n]+)/ms;
	($statusHash->{nodes}) 		= $jobcheck =~ /Allocated Nodes:\s*\n\s*(\S+?)/ms;
	($statusHash->{nodecount}) 	= $jobcheck =~ /NodeCount:\s*(\d+)/ms;


	($statusHash->{status}) 	||= '';
	($statusHash->{submittime}) ||= '';
	($statusHash->{starttime})  ||= '';
	($statusHash->{nodes})  	||= '';
	($statusHash->{nodecount})  ||= '';


	return $statusHash;	
}


=head

    SUBROUTINE      statusHash    

    PURPOSE

		HASH JOBID AGAINST STATUS

=cut

method statusHash {


	my $joblines = $self->jobLines();
	return {} if not defined $joblines;

	my $statusHash;
	foreach my $line ( @$joblines )
	{

		if ( $line =~ /^\s*(\d+)/ )
		{
			$statusHash->{$1} = $self->lineStatus($line);
		}
	}


	return $statusHash;
}


=head2

	SUBROUTINE		lineStatus

	PURPOSE

		RETURN THE STATUS ENTRY FOR A JOB LINE

	NOTES

		SGE JOB LIFECYCLE:

		"qw"	queued and waiting
		"t"		transferring to an available node
		"r" 	running
		"Eqw"	error state

		When a job no longer appears on the qstat output, it has finished or been deleted.

		QSTAT FORMAT:
		job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID 
		-----------------------------------------------------------------------------------------------------------------
			135 0.50000 tophatBatc www-data     Eqw   04/13/2011 00:39:46                                    3 1
			136 0.50000 tophatBatc www-data     Eqw   04/13/2011 01:14:03                                    3 1
			170 0.00000 tophatBatc root         qw    04/13/2011 17:00:25                                    3 1      

=cut

method lineStatus ($line) {	
	my ($status) = $line =~ /^\s*\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/;

	return "running" if $status eq "r";
	return "starting" if $status eq "t";
	return "queued" if $status eq "qw";
	return "error" if $status eq "Eqw";
}




} #### Agua::Monitor::SGE



1;




