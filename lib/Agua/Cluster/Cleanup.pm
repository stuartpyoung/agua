package Agua::Cluster::Cleanup;
use Moose::Role;
use Moose::Util::TypeConstraints;

use Data::Dumper;


################################################################################
##########################      CLEANUP METHODS       ##########################
################################################################################
=head2

	SUBROUTINE     deleteMiscfiles

    PURPOSE

        REMOVE THE ALIGNMENT SUBDIRS USED TO PERFORM INPUT FILE CHUNK ALIGNMENTS

    INPUT

        1. LOCATION OF OUTPUTDIR 

=cut

sub deleteMiscfiles {
	my $self		=	shift;
	my $directory	=	shift;
	my $references	=	shift;

	#### CHECK FOR OUTPUT DIRECTORY
	print "Agua::Cluster::Cleanup::deleteMiscfiles    directory not defined: $directory\n" if not defined $directory or not $directory;
	print "Agua::Cluster::Cleanup::deleteMiscfiles    directory is '': $directory\n" if not defined $directory;
	print "Agua::Cluster::Cleanup::deleteMiscfiles    Can't find directory: $directory\n" if not -d $directory;

	##### CHDIR TO OUTPUT DIRECTORY
	#chdir($directory) or die "Could not CHDIR to directory: $directory\n";
	#opendir(DIR, $directory) or die "Can't open directory: $directory\n";
	#my @subdirs = readdir(DIR);
	#close(DIR);
	#
	#foreach my $subdir ( @subdirs )
	#{
	#	next if not $subdir =~ /^\d+$/;
	#	
	#	foreach my $reference ( @$references )
	#	{
	#		my $refdir = "$directory/$subdir/$reference";
	#		opendir(REFDIR, $refdir) or die "Can't open refdir: $refdir\n";
	#		my @reffiles = readdir(REFDIR);
	#		close(REFDIR);
	#		
	#		foreach my $reffile ( @reffiles )
	#		{
	#			next if not $reffile =~ /^\d+$/;
	#			`rm -fr $refdir/$reffile`;
	#		}
	#	}
	#}

}



=head2

	SUBROUTINE     archiveMiscfiles

    PURPOSE

        REMOVE THE ALIGNMENT SUBDIRS USED TO PERFORM INPUT FILE CHUNK ALIGNMENTS

    INPUT

        1. LOCATION OF OUTPUTDIR 

=cut

sub archiveMiscfiles {
	my $self		=	shift;
	my $directory	=	shift;
	my $references	=	shift;

	#### CHECK FOR OUTPUT DIRECTORY
	print "Agua::Cluster::Cleanup::archiveMiscfiles    directory not defined: $directory\n" if not defined $directory or not $directory;
	print "Agua::Cluster::Cleanup::archiveMiscfiles    directory is '': $directory\n" if not defined $directory;
	print "Agua::Cluster::Cleanup::archiveMiscfiles    Can't find directory: $directory\n" if not -d $directory;

	##### CHDIR TO OUTPUT DIRECTORY
	#chdir($directory) or die "Could not CHDIR to directory: $directory\n";
	#opendir(DIR, $directory) or die "Can't open directory: $directory\n";
	#my @subdirs = readdir(DIR);
	#close(DIR);
	#
	#foreach my $subdir ( @subdirs )
	#{
	#	next if not $subdir =~ /^\d+$/;
	#	
	#	foreach my $reference ( @$references )
	#	{
	#		my $refdir = "$directory/$subdir/$reference";
	#		opendir(REFDIR, $refdir) or die "Can't open refdir: $refdir\n";
	#		my @reffiles = readdir(REFDIR);
	#		close(REFDIR);
	#		
	#		foreach my $reffile ( @reffiles )
	#		{
	#			next if not $reffile =~ /^\d+$/;
	#			`rm -fr $refdir/$reffile`;
	#		}
	#	}
	#}	
}

=head2

	SUBROUTINE     deleteAlignmentSubdirs

    PURPOSE

        REMOVE THE ALIGNMENT SUBDIRS USED TO PERFORM INPUT FILE CHUNK ALIGNMENTS

    INPUT

        1. LOCATION OF OUTPUTDIR 

=cut

sub deleteAlignmentSubdirs {
	my $self		=	shift;
	my $directory	=	shift;
	my $references	=	shift;

	#### CHECK FOR OUTPUT DIRECTORY
	print "Agua::Cluster::Cleanup::deleteAlignmentSubdirs    directory not defined: $directory\n" if not defined $directory or not $directory;
	print "Agua::Cluster::Cleanup::deleteAlignmentSubdirs    directory is '': $directory\n" if not defined $directory;
	print "Agua::Cluster::Cleanup::deleteAlignmentSubdirs    Can't find directory: $directory\n" if not -d $directory;

	#### CHDIR TO OUTPUT DIRECTORY
	chdir($directory) or die "Could not CHDIR to directory: $directory\n";
	opendir(DIR, $directory) or die "Can't open directory: $directory\n";
	my @subdirs = readdir(DIR);
	close(DIR);
	print "subdirs: @subdirs\n";

	foreach my $subdir ( @subdirs )
	{
		next if not $subdir =~ /^\d+$/;
		print "subdir: $subdir\n";

		foreach my $reference ( @$references )
		{
			my $refdir = "$directory/$subdir/$reference";
			opendir(REFDIR, $refdir) or die "Can't open refdir: $refdir\n";
			my @reffiles = readdir(REFDIR);
			close(REFDIR);

			foreach my $reffile ( @reffiles )
			{
				next if not $reffile =~ /^\d+$/;
				print "rm -fr $refdir/$reffile\n";
				`rm -fr $refdir/$reffile`;
			}
		}		
	}
}

=head2

	SUBROUTINE     deleteSplitfiles

    PURPOSE

        REMOVE THE INPUT FILE CHUNKS (SPLITFILES)

    INPUT

        1. LOCATION OF OUTPUTDIR 

=cut

sub deleteSplitfiles {
	my $self		=	shift;
	my $directory	=	shift;

	#### CHECK FOR OUTPUT DIRECTORY
	print "Agua::Cluster::Cleanup::deleteSplitfiles    directory not defined: $directory\n" and exit if not defined $directory or not $directory;
	print "Agua::Cluster::Cleanup::deleteSplitfiles    Can't find directory: $directory\n" and exit if not -d $directory;

	#### CHDIR TO OUTPUT DIRECTORY
	chdir($directory) or die "Could not CHDIR to directory: $directory\n";
	opendir(DIR, $directory) or die "Can't open directory: $directory\n";
	my @subdirs = readdir(DIR);
	close(DIR);
	foreach my $subdir ( @subdirs )
	{
		next if not $subdir =~ /^\d+$/;
		opendir(SUBDIR,"$directory/$subdir") or die "Can't open subdir: $directory/$subdir\n";
		my @subfiles = readdir(SUBDIR);
		close(SUBDIR);
		foreach my $subfile ( @subfiles )
		{
			next if not $subfile =~ /^\d+$/;
			print "subfile: $directory/$subdir/$subfile\n";
			`rm -fr $directory/$subdir/$subfile`;
		}
	}
}




1;
