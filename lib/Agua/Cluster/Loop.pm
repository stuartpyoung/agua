package Agua::Cluster::Loop;
use Moose::Role;


=head2

	PACKAGE		Agua::Cluster::Loop

	PURPOSE

		1. REPEATEDLY EXECUTE AN APPLICATION, CHANGING THE VALUE OF

			THE SPECIFIED PARAMETER EVERY TIME AND FOR A SPECIFIED NUMBER

			OF REPLICATES

		2. DO THE SAME FOR AN ARRAY OF WORKFLOW APPLICATIONS, CHANGING

			THE VALUES FOR DOWNSTREAM APPLICATIONS EACH TIME

	VERSION

		0.02	Added '--loop' option with args 'serial', 'parallel' and 'distributed'

		O.01	Basic loop with serial or 'concurrent'

=cut

use strict;
use warnings;

# STRINGS
has 'arguments'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'executable'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'parameter'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'values'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'replicates'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'concurrent'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'loopout'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'label'			=> ( isa => 'Str|Undef', is => 'rw', default => '' );

#### EXTERNAL MODULES
use Data::Dumper;
use File::Path;
use threads qw(yield);

=head2

	SUBROUTINE		runLoop

	PURPOSE

		RUN THE APPLICATION FOR THE DESIGNATED PERMUTATION OF LOOPS

=cut

sub runLoop {
	my $self		=	shift;

	#### GET ARGUMENT VALUES
	my $executable 	=	$self->loop_value("executable");
	my $parameter 	= 	$self->loop_value("parameter");
	my $values 		= 	$self->loop_value("values");
	my $replicates 	= 	$self->loop_value("replicates");
	my $loop 		= 	$self->loop_value("loop");
	my $loopout 	= 	$self->loop_value("loopout");
	my $label 		=	$self->loop_value("looplabel");

	#### EXTRACT CLUSTER-RELATED ARGUMENTS IF loop IS distributed
	my ($cluster, $outdir, $maxjobs, $queue);
	if ( $loop eq "distributed" )
	{
		$cluster=	$self->loop_value("cluster");
		$outdir =	$self->loop_value("outdir");
		$maxjobs=	$self->loop_value("maxjobs");
		$queue	=	$self->loop_value("queue");

		$self->loop_value("walltime");
		$self->loop_value("cpus");
		$self->loop_value("qstat");
		$self->loop_value("qsub");
		$self->loop_value("sleep");
		$self->loop_value("cleanup");
		$self->loop_value("dot");
	}

	my $oldout;
	if ( defined $loopout )
	{
		print "Agua::Logic::Loop::runLoop    Printing STDOUT to loopout:\n\n$loopout\n\n";

		open $oldout, ">&STDOUT" or die "Can't open old STDOUT\n";
		#open $olderr, ">&STDERR" or die "Can't open old STDERR\n";

		my ($loopout_path) = $loopout =~ /^(.+)\/[^\/]+$/;
		File::Path::mkpath($loopout_path) if not -d $loopout_path;
		open(STDOUT, ">$loopout") or die "Can't open STDOUT file: $loopout\n" if defined $loopout;
	}

	#### GET THE REMAINING ARGUMENTS FOR THE EXECUTABLE
	my $arguments	=	$self->get_arguments();

	#### DEBUG
	print "Agua::Logic::Loop::runLoop     executable: $executable\n";
	print "Agua::Logic::Loop::runLoop     parameter:  $parameter\n" if defined $parameter;
	print "Agua::Logic::Loop::runLoop     values:     $values\n" if defined $values;
	print "Agua::Logic::Loop::runLoop     replicates: $replicates\n" if defined $replicates;
	print "Agua::Logic::Loop::runLoop     loop: $loop\n" if defined $loop;
	print "Agua::Logic::Loop::runLoop     cluster: $cluster\n" if defined $cluster;
	print "Agua::Logic::Loop::runLoop     outdir: $outdir\n" if defined $outdir;
	print "Agua::Logic::Loop::runLoop     label: $label\n" if defined $label;

	#### CHECK loop IS DEFINED AND CORRECT VALUE
	print "Agua::Logic::Loop::runLoop    loop not defined. Exiting\n" and exit if not defined $loop;
	if ( $loop eq "distributed" )
	{
		print "Agua::Logic::Loop::runLoop    cluster not defined. Exiting\n" and exit if not defined $cluster;
		print "Agua::Logic::Loop::runLoop    label not defined. Exiting\n" and exit if not defined $label;
		print "Agua::Logic::Loop::runLoop    outdir not defined. Exiting\n" and exit if not defined $outdir;
		print "Agua::Logic::Loop::runLoop    maxjobs not defined. Exiting\n" and exit if not defined $maxjobs;
		print "Agua::Logic::Loop::runLoop    queue not defined. Exiting\n" and exit if not defined $queue;
	}

	print "Agua::Logic::Loop::runLoop    loop value not supported: $loop\n" and exit
		if $loop !~ /^(serial|parallel|distributed)$/;

	print "Agua::Logic::Loop::runLoop    Neither parameter nor replicates is defined. Exiting.\n" if not defined $parameter and not defined $replicates;
	print "Agua::Logic::Loop::runLoop    values not defined for parameter. Exiting.\n" if defined $parameter and not defined $values;

	#### GET VALUES AND REPLICATES ARRAYS
	my $values_array = [];
	@$values_array = split ",", $values if defined $values;

	my $replicates_array = $self->stringToArray($replicates);

	#### GET LOOPED COMMANDS
	my $commands = $self->loop($executable, $arguments, $parameter, $values_array, $replicates_array);
#exit;

	#### COLLECT OVERALL STATUS AND FAILED OR DUBIOUS JOB LABELS
	my $status = 'unknown';
	my $sublabels = '';

	#### RUN COMMANDS ON CLUSTER
	if ( $loop eq "distributed" )
	{
		print "Agua::Logic::Loop::runLoop    Running loop in $loop mode\n";

		my $jobs = [];
		for ( my $i = 0; $i < @$commands; $i++ )
		{
			my $replicate_number = $$replicates_array[$i];

			my $command = $$commands[$i];
			my $joblabel = "$label-$replicate_number";
			my $scriptfile = "$outdir/$label-$replicate_number.sh";
			my $usagefile = "$outdir/$label-$replicate_number-usage.txt";
			my $stdoutfile = "$outdir/$label-$replicate_number-stdout.txt";
			my $stderrfile = "$outdir/$label-$replicate_number-stderr.txt";

			my $job = $self->setJob([$command], $joblabel, $outdir, $scriptfile, $usagefile, $stdoutfile, $stderrfile);
			push @$jobs, $job;
		}

		#### USE LIB FOR CLUSTER MONITOR
		print "Agua::Logic::Loop::runLoop    cluster: **$cluster**\n";

		if ( $cluster =~ /^PBS$/ )
		{
			print "Agua::Logic::Loop::runLoop    DOING require Monitor::PBS\n";
			eval "require Monitor::PBS";
		}
		elsif ( $cluster =~ /^LSF$/ )
		{
			print "Agua::Logic::Loop::runLoop    DOING require Monitor::LSF\n";
			eval "require Monitor::LSF";
		}
		elsif ( $cluster =~ /^SGE$/ )
		{
			print "Agua::Logic::Loop::runLoop    DOING require Monitor::LSF\n";
			eval "require Monitor::SGE";
		}
		else
		{
			print "Agua::Logic::Loop::runLoop    cluster $cluster did not match LSF or PBS\n";
		}

		#### WILL AUTOMATICALLY CHECK OUTPUTS FOR COMPLETION
		#### AND PRINT STATUS SIGNAL TO STDOUT 
		($status, $sublabels) = $self->runJobs($jobs, $label);

		#### PRINT JOB STATUS SIGNAL
		print "\n------------------------------------------------------------\n";
		print "---[status $label: $status $sublabels]---";
		print "\n------------------------------------------------------------\n";		

		print "Agua::Logic::Loop::runLoop    Finished doing distributed jobs\n";
	}

	#### RUN COMMANDS IN PARALLEL LOCALLY
	elsif ( $loop eq "parallel" )
	{
		print "Agua::Logic::Loop::runLoop    DOING CONCURRENT JOBS\n";
		my $threads = [];
		for my $command ( @$commands )
		{
			print "Agua::Logic::Loop::runLoop    CONCURRENT: $command\n";
			my $thread = threads->new(
				sub {
					return `$command`;
				}
			);
			sleep(1);
			push @$threads, $thread;
		}

		my $outputs = [];
		foreach my $thread ( @$threads )
		{
			my $output = $thread->join;
			print "Agua::Logic::Loop::runLoop    OUTPUT: $output\n";
			push @$outputs, $output;
		}

		#### CHECK OUTPUTS FOR COMPLETION
		($status, $sublabels) = $self->completionStatus($outputs);

		#### PRINT JOB STATUS SIGNAL
		print "\n------------------------------------------------------------\n";
		print "---[status $label: $status $sublabels]---";
		print "\n------------------------------------------------------------\n";		

		print "Agua::Logic::Loop::runLoop    Finished doing parallel jobs\n";
	}

	#### RUN COMMANDS IN SERIES
	elsif ( $loop eq "serial" )
	{
		print "Agua::Logic::Loop::runLoop    Running loop in 'serial' mode\n";

		#### RUN THE COMMANDS
		my $outputs = [];
		foreach my $command ( @$commands )
		{
			print "Agua::Logic::Loop::runLoop     $command\n";
			my $output = `$command`;
			push @$outputs, $output;
		}

		#### CHECK OUTPUTS FOR COMPLETION 
		($status, $sublabels) = $self->completionStatus($outputs);

		($label) = $executable =~ /([^\/]+)$/ if not defined $label;

		#### PRINT JOB STATUS SIGNAL
		print "\n------------------------------------------------------------\n";
		print "---[status $label: $status $sublabels]---";
		print "\n------------------------------------------------------------\n";		

		print "Agua::Logic::Loop::runLoop    Finished doing distributed jobs\n";
	}
	else
	{
		print "Agua::Logic::Loop::runLoop    Exiting. Missing handler for loop mode: $loop\n";
		exit;
	}


	##### RESTORE OLD STDOUT 
	if ( defined $oldout )
	{
		print "Agua::Logic::Loop::runLoop    Redirecting STDOUT back to standard output\n";
		open STDOUT, ">&", $oldout;
	}

	#### PRINT JOB STATUS SIGNAL
	print "\n------------------------------------------------------------\n";
	print "---[status $label: $status $sublabels]---";
	print "\n------------------------------------------------------------\n";
}

sub completionStatus {
	my $self		=	shift;
	my $outputs		=	shift;

	$outputs = [$outputs] if ref($outputs) ne "ARRAY";

	#### COLLECT JOB COMPLETION SIGNAL (AND SUBLABELS OF INCOMPLETE
	#### JOBS, MISSING FILES OR BOTH - I.E., FAILED JOBS)
	#### STATUS REPORTING HIERARCHY: completed < incomplete < missing < failed
	my $overall_status = "completed";
	my $overall_sublabels = '';
	my ($label, $status, $sublabels);
	foreach my $output (@$outputs)
	{
		if ( $output =~ /---\[status\s+(\S+):\s+(\S+)\s*(\S*)\]/ms )
		{
			$label = $1;
			$status = $2;
			$sublabels = $3;
			print "Job label '$label' completion signal: $status\n";
			$overall_status = "complete" if $status eq "complete"
				and $overall_status ne "incomplete"
				and $overall_status ne "missing"
				and $overall_status ne "failed";
			$overall_status = "incomplete" if $status eq "incomplete"
				and $overall_status ne "missing"
				and $overall_status ne "failed";
			$overall_status = "missing" if $overall_status eq "missing"
				and $overall_status ne "failed";
			$overall_status = "failed" if $status eq "failed";
			$overall_sublabels .= $sublabels . "," if $sublabels
				and $status ne "complete";
		}
	}

	return ($overall_status, $overall_sublabels);
}



=head2

	SUBROUTINE		loop

	PURPOSE

        1. REPEATEDLY EXECUTE AN APPLICATION, CHANGING THE VALUE OF

			THE SPECIFIED PARAMETER EVERY TIME

		2. REPEAT FOR A SPECIFIED NUMBER OF REPLICATES

    INPUT

        1. EXECUTABLE AND ITS ARGUMENTS (WITH STRING '%VALUE%' 

			IN THE PLACE OF THE ACTUAL PARAMETER VALUE)

		2. PARAMETER TO BE CHANGED

		3. COMMA-SEPARATED LIST OF VALUES FOR THE PARAMETER

		4. COMMA-SEPARATED LIST OF REPLICATES

    OUTPUT

        1. OUTPUTS OF EACH RUN OF THE EXECUTABLE USING A

			DIFFERENT VALUE FOR THE PARAMETER EACH TIME

=cut

sub loop {

	my $self		=	shift;
	my $executable	=	shift;
	my $arguments	=	shift;
	my $parameter	=	shift;
	my $values		=	shift;
	my $replicates	=	shift;

	print "Agua::Logic::Loop::loop     executable: $executable\n";
	print "Agua::Logic::Loop::loop     arguments:  @$arguments\n";
	print "Agua::Logic::Loop::loop     parameter:  $parameter\n" if defined $parameter;
	print "Agua::Logic::Loop::loop     values:     @$values\n" if defined $values;
	print "Agua::Logic::Loop::loop     replicates: @$replicates\n" if defined $replicates;

	my $commands;

	#### RUN replicates TIMES WITH values VALUES
	for ( my $counter = 0; $counter < @$replicates; $counter++ )
	{
		my $replicate = $$replicates[$counter];

		if ( defined $parameter )
		{
			for ( my $i = 0; $i < @$values; $i++ )
			{
				my $instance_args;
				@$instance_args = @$arguments;

				my $value = $$values[$i];

				#### SUBSTITUTE replicate FOR ONE OR MORE '%REPLICATE%' STRINGS IN ALL ARGUMENTS
				$instance_args = $self->fill_in($instance_args, "%REPLICATE%", $replicate);

				#### SUBSTITUTE value FOR ONE OR MORE '%VALUE%' STRINGS IN ALL ARGUMENTS
				$instance_args = $self->fill_in($instance_args, "%VALUE%", $value);

				#### SUBSTITUTE parameter FOR ONE OR MORE '%PARAMETER%' STRINGS IN ALL ARGUMENTS
				$instance_args = $self->fill_in($instance_args, "%PARAMETER%", $parameter);

				#### ADD parameter ARGUMENT TO FRONT OF ARGS
				unshift @$instance_args, "$parameter $value";

				my $command = "$executable @$instance_args";

				push @$commands, $command;
			}
		}
		else
		{
			my $instance_args;
			@$instance_args = @$arguments;

			#### SUBSTITUTE replicate FOR ONE OR MORE '%REPLICATE%' STRINGS IN ALL ARGUMENTS
			$instance_args = $self->fill_in($instance_args, "%REPLICATE%", $replicate);

			foreach my $arg ( @$instance_args )
			{
				$arg = qq{"$arg"};
			}

			my $command = "$executable @$instance_args";

#exit;

			push @$commands, $command;
		}
	}

	return $commands;
}

=head2

	SUBROUTINE		set_parameter

	PURPOSE

		SET THE VALUE OF A PARAMETER IN arguments

=cut
sub set_parameter {
	my $self		=	shift;
	my $arguments		=	shift;
	my $parameter			=	shift;
	my $value			=	shift;


	for ( my $i = 0; $i < @$arguments; $i++ )
	{
		if ( "--$parameter" eq $$arguments[$i] )
		{
			$$arguments[$i + 1] = $value;
			return $arguments;
		}	
	}

	return $arguments;
}



=head2

	SUBROUTINE		fill_in

	PURPOSE

		SUBSTITUTE counter FOR ONE OR MORE '%REPLICATE%' STRINGS IN ALL ARGUMENTS

=cut
sub fill_in {
	my $self		=	shift;
	my $arguments		=	shift;
	my $pattern			=	shift;
	my $value			=	shift;

	foreach my $argument ( @$arguments )
	{
		$argument =~ s/$pattern/$value/ig;
	}

	return $arguments;
}


=head2

	SUBROUTINE		loop_value

	PURPOSE

		1. EXTRACT THE VALUE FOR THE NAMED ARGUMENT FROM THE arguments

			ARRAY AND RETURN IT

		2. REPLACE THE EXITING arguments ARRAY WITH THE NEWLY TRUNCATED ONE

=cut

sub loop_value {
	my $self			=	shift;
	my $name			=	shift;


	my $arguments		=	$self->get_arguments();

	my $value;
	for ( my $i = 0; $i < @$arguments; $i++ )
	{
		if ( "--$name" eq $$arguments[$i] )
		{
			$value = $$arguments[$i + 1];
			splice @$arguments, $i, 2;
			$self->{"_$name"} = $value;
		}
	}


	$self->set_arguments($arguments);
	#$self->{_arguments} = $arguments;

	return $value;
}







=head2

	SUBROUTINE		oldLoop

	PURPOSE

        1. REPEATEDLY EXECUTE AN APPLICATION, CHANGING THE VALUE OF

			A PARTICULAR PARAMETER EVERY TIME

		2. REPEAT FOR A SPECIFIED NUMBER OF REPLICATES

    INPUT

        1. EXECUTABLE AND ITS ARGUMENTS (WITH STRING '%VALUE%' 

			IN THE PLACE OF THE ACTUAL PARAMETER VALUE)

		2. PARAMETER TO BE CHANGED

		3. COMMA-SEPARATED LIST OF VALUES FOR THE PARAMETER

		4. COMMA-SEPARATED LIST OF REPLICATES

    OUTPUT

        1. OUTPUTS OF EACH RUN OF THE EXECUTABLE USING A

			DIFFERENT VALUE FOR THE PARAMETER EACH TIME

=cut

sub oldLoop {
	my $self		=	shift;
	my $executable	=	shift;
	my $arguments	=	shift;
	my $parameter	=	shift;
	my $values		=	shift;
	my $replicates	=	shift;

	#### RUN replicates TIMES WITH values VALUES
	for ( my $counter = 0; $counter < @$replicates; $counter++ )
	{
		my $replicate = $$replicates[$counter];

		for ( my $i = 0; $i < @$values; $i++ )
		{
			my $instance_args;
			@$instance_args = @$arguments;

			my $value = $$values[$i];

			#### SUBSTITUTE replicate FOR ONE OR MORE '%REPLICATE%' STRINGS IN ALL ARGUMENTS
			$instance_args = fill_in($instance_args, "%REPLICATE%", $replicate);

			#### SUBSTITUTE value FOR ONE OR MORE '%VALUE%' STRINGS IN ALL ARGUMENTS
			$instance_args = fill_in($instance_args, "%VALUE%", $value);

			#### SUBSTITUTE parameter FOR ONE OR MORE '%PARAMETER%' STRINGS IN ALL ARGUMENTS
			$instance_args = fill_in($instance_args, "%PARAMETER%", $parameter);

			my $command = "$executable @$instance_args";
			print "command: $command\n\n\n";

			`$command`;		
		}
	}
}



1;

