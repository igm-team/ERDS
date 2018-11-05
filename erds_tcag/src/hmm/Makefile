#
#		Authors: Mingfu Zhu & David B. Goldstein
#		Organization: Center for Human Genome Variation, Duke School of Medicine
#		Modified from Tapas Kanungo, "UMDHMM: Hidden Markov Model Toolkit," in "Extended Finite
#		State Models of Language," A. Kornai (editor), Cambridge University Press, 1999.																	    
#

CFLAGS= -g
INCS=
# use the following line to "Purify" the code
#CC=purify gcc
CC=gcc
SRCS=baum.c viterbi.c forward.c backward.c hmmutils.c nrutil.c

all :hmm

hmm: hmm.o backward.o baum.o forward.o hmmutils.o nrutil.o viterbi.o
	 $(CC) -o hmm hmm.o backward.o baum.o forward.o hmmutils.o nrutil.o viterbi.o -lm
clean:
	rm *.o a.out 
# DO NOT DELETE THIS LINE -- make depend depends on it.

