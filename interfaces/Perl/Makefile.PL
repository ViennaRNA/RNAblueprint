# File : Makefile.PL
use ExtUtils::MakeMaker;
WriteMakefile(
	`NAME'    => `example',                  # Name of package
	`LIBS'    => [`-lm'],                    # Name of custom libraries
	`OBJECT'  => `example.o example_wrap.o'  # Object files
);


