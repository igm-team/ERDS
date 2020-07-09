#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <unistd.h>
#include <sys/stat.h>

using namespace std;
const int MAX_LINE=128*1024*1024;
const int BAD_COMMAND_LINE_ARGS_FLAG = -1;

vector<int>* rd_bp = new vector<int>();
vector<int>* rcv = new vector<int>();
vector<int>* mqv = new vector<int>();

// From http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
bool file_exists (const string& name) {
	struct stat buffer;   
	return (stat (name.c_str(), &buffer) == 0); 
}

// The following two functions are from http://stackoverflow.com/questions/236129/split-a-string-in-c
// They emulate the functionality of the perl "split" function
template<typename Out>
void split(const string &s, char delim, Out result) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        *(result++) = item;
    }
}

vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, back_inserter(elems));
    return elems;
}


double midhalf(vector<int>* v) { // sub midhalf{
	int count = v->size(); // my $count = $#array+1;
	long int sum = 0; // my $sum;
	vector<int>* array = v;
	sort(array->begin(), array->end()); // my @sarray = sort { $a <=> $b } @array;
	int upper = min(count - 1, (int) (count * 0.75 + 0.5));
	int bottom = count * 0.25 + 0.5;
	for (int i = bottom; i <= upper; i++) { // for(my $i=$bottom; $i<=$upper; $i++){
		sum += array->at(i); // $sum+=$sarray[$i];
	}
	return (double) (sum * count / (upper - bottom + 1)); // return $sum*$count/($upper-$bottom+1);
}

int count_CIGAR(string CIGAR_str, bool use_I) {
	int start_pos = -1;
	int total = 0;
	for (int i = 0; i < CIGAR_str.length(); i++) {
		int ascii = (int) CIGAR_str.at(i);
		
		if ((ascii == 77 || (ascii == 73 && use_I)) && i != 0 && start_pos != -1) { // 77 -> M, 73 -> I
			int len = i - start_pos;
			string to_convert = CIGAR_str.substr(start_pos, len);
			total += stoi(to_convert);
			start_pos = -1;
		} else if (ascii >= 48 && ascii <= 57) {
			if (start_pos == -1) {
				start_pos = i;
			}
		} else {
			start_pos = -1;
		}
	}
	
	return total;
}

// sub rc{
// my ($start, $alg_length, $mq)=@_;
void rc(int start, int alg_length, int mq, int win_size_del) {
	int end = start + alg_length - 1; // my $end=$start+$alg_length-1;
	for (size_t i = start; i <= end; i++) { // for(my $i=$start; $i<=$end; $i++){
		// Make sure rd_bp[i] actually exists; if not, create it
		if (rd_bp->size() < i + 1) {
			rd_bp->resize(i + 1);
			rd_bp->at(i) = 0;
		}     
		rd_bp->at(i) = rd_bp->at(i) + 1;
	}

	int s = (start - 1) / win_size_del; // my $s=int(($start-1)/$win_size_del);
	int e = (end - 1) / win_size_del; // my $e=int(($end-1)/$win_size_del);
	// Increase size of rcv and mqv if needed.
	int old_size = rcv->size();
	if (old_size < e + 1) {
		rcv->resize(e + 1, -1);
		mqv->resize(e + 1, -1);
	}
	
	#ifdef DEBUG_READ_INFO
		cout << "s: " << s << "\n";
		cout << "e: " << e << "\n";
	#endif

	if (rcv->at(e) == -1) {
		rcv->at(e) = 0;
	}

	if (mqv->at(e) == -1) {
		mqv->at(e) = 0;
	}

	if (s == e) { // if($s==$e) {
		rcv->at(s) = rcv->at(s) + alg_length; // $rc[$s]+=$alg_length;
		mqv->at(s) = mqv->at(s) + (mq * alg_length); // $mq[$s]+=$mq*$alg_length;
	} else {
		
		if (rcv->at(s) == -1) {
			rcv->at(s) = 0;
		}

		if (mqv->at(s) == -1) {
			mqv->at(s) = 0;
		}

		rcv->at(s) = rcv->at(s) + (win_size_del * (s + 1) - start + 1); 		// $rc[$s]+=($win_size_del*($s+1)-$start+1);
		rcv->at(e) = rcv->at(e) + (end - win_size_del * e);			// $rc[$e]+=($end-$win_size_del*$e);
		mqv->at(s) = mqv->at(s) + (mq * (win_size_del * (s + 1) - start + 1));	// $mq[$s]+=$mq*($win_size_del*($s+1)-$start+1);
		mqv->at(e) = mqv->at(e) + (mq * (end - win_size_del * e));		// $mq[$e]+=$mq*($end-$win_size_del*$e);
		
		for (int i = s + 1; i <= e - 1; i++) {			// for(my $i=$s+1; $i<=$e-1; $i++){
			if (rcv->at(i) == -1) {
				rcv->at(i) = 0;
			}
			if (mqv->at(i) == -1) {
				mqv->at(i) = 0;
			}
			rcv->at(i) = rcv->at(i) + win_size_del;				// $rc[$i]+=$win_size_del;
			mqv->at(i) = mqv->at(i) + (mq * win_size_del);			// $mq[$i]+=$mq*$win_size_del;
		}
	}
}

int main(int argc, char **argv) {

	if (argc != 13) {
		cout << "Input error: incorrect number of command-line arguments.\n";
		cout << "Usage: " << argv[0] << " bam_file ins_pem min_det_size chr aligned_length_thh CIGAR_length large_pem win_size_del read_length output_file\n";
		return BAD_COMMAND_LINE_ARGS_FLAG;
	}

	string bam_file 	= argv[1];       // $bam_file=$_[0];
	int ins_pem 		= stoi(argv[2]); // Global variable in perl version; calculated in test_bam
	int min_det_size 	= stoi(argv[3]); // Global variable in perl version; calculated in test_bam
	string chr 		= argv[4];       // Global variable in perl version
	int aligned_length_thh 	= stoi(argv[5]); // Global variable (from parameter file) in perl version; default = 30
	int CIGAR_length 	= stoi(argv[6]); // Global variable (from parameter file) in perl version; default = 12
	int large_pem		= stoi(argv[7]); // Global variable (from parameter file) in perl version; default = 50000
	int win_size_del	= stoi(argv[8]); // Global variable (from parameter file) in perl version; default = 200
	int read_length		= stoi(argv[9]); // Global variable in perl version; calculated in test_bam
	string chr_rc		= argv[10];
	string chr_sd_rc	= argv[11];
	string output_file	= argv[12];      // Not needed in Perl version because file is opened and closed outside of read_bam function

	setvbuf(stdout, NULL, _IONBF, 0); // Stop buffering stdout, so even if the program crashes, we will get to see any output up to that point
	
	int ins_cut = ins_pem + min_det_size; // my $ins_cut=$ins_pem+$min_det_size;

	cout << "### PARAMETERS ###\n";
	cout << "bam_file: " << bam_file << "\n";
	cout << "ins_pem: " << ins_pem << "\n";
	cout << "min_det_size: " << min_det_size << "\n";
	cout << "chr: " << chr << "\n";
	cout << "aligned_length_thh: " << aligned_length_thh << "\n";
	cout << "CIGAR_length: " << CIGAR_length << "\n";
	cout << "large_pem: " << large_pem << "\n";
	cout << "win_size_del: " << win_size_del << "\n";
	cout << "read_length: " << read_length << "\n";
	cout << "chr_rc: " << chr_rc << "\n";
	cout << "chr_sd_rc: " << chr_sd_rc << "\n";
	cout << "output_file: " << output_file << "\n";

	string touch_command = "touch '" + chr_sd_rc + "'";
	system(touch_command.c_str());
	
	vector<string>* rowv = new vector<string>();
	vector<int>* posv = new vector<int>();
	vector<int>* pemv = new vector<int>(); // my (@row,@pos,@pem);

	string samtools_command = "samtools view -F 2048 '" + bam_file + "' " + chr;
	cout << "samtools command is [" << samtools_command << "]\n"; // FOR DEBUGGING

	ofstream PSAM;
	PSAM.open(output_file, ios::out | ios::app); // open(PSAM,">>$chr_pem_sam") or die "Can't open $chr_pem_sam";

	FILE *sam_stream;
	char *buff = new char[MAX_LINE];

	if (!(sam_stream = popen(samtools_command.c_str(), "r"))) {
	        return -1;
	}
	
	string in;

	while (!feof(sam_stream)) {
		if (fgets(buff, MAX_LINE, sam_stream) != NULL) { // while(my $in=<ST>) {
				
			in = buff;
			in.erase(remove(in.begin(), in.end(), '\n'), in.end()); // $in=~s/[\r\n]//g; // from http://stackoverflow.com/questions/1488775/c-std::remove-new-line-from-multiline-string
			in.erase(remove(in.begin(), in.end(), '\r'), in.end()); // $in=~s/[\r\n]//g; // from http://stackoverflow.com/questions/1488775/c-std::remove-new-line-from-multiline-string
			vector<string>* list = new vector<string>();
			*list = split(in, '\t'); // @list=split /\s/, $in;
			int pos, mq;

			try {
				pos = stoi(list->at(3)); // $pos=$list[3];
				mq  = stoi(list->at(4)); // $mq=$list[4];
			} catch (const std::exception & e){
				cerr << "Possible long line" << endl;
			        continue;
			}

			string CIGAR = list->at(5); // $CIGAR=$list[5];
			int alg_length = 0; // $alg_length=0;
			
			// while($CIAGR=~/(\d+)[M]/g){
				// $alg_length+=$1;
			// }
			alg_length = count_CIGAR(CIGAR, false);

			// next if($alg_length<$aligned_length_thh or $CIAGR eq "*" or length($CIAGR)>$CIAGR_length);
			if (alg_length < aligned_length_thh || CIGAR.compare("*") == 0 || CIGAR.length() > CIGAR_length) {
				continue;
			}

			int ins_read = abs(stoi(list->at(8))); // $ins_read=abs($list[8]);
			string bq = list->at(10); // $bq=$list[10];
			list->at(5) = "1M"; // $list[5]="1M";
			
			string row = "";

			if (ins_read > ins_cut && ins_read < 2*large_pem) { // if($ins_read>$ins_cut and $ins_read<2*$large_pem){
				
				// $row="$list[0]";
				// for(my $j=1; $j<=8; $j++){
				// 	$row=$row."\t$list[$j]";
				// }
				// $row=$row."\tA\tH\tMD:Z:$CIAGR";
				row = list->at(0) + "\t" + list->at(1) + "\t" + list->at(2) + "\t" + list->at(3) + "\t" + list->at(4) + "\t" + list->at(5) + "\t" + list->at(6) + "\t" + list->at(7) + "\t" + list->at(8) + "\tA\tH\tMD:Z:" + CIGAR;
				rowv->push_back(row); // push(@row,$row);
				posv->push_back(pos); // push(@pos,$pos);
				pemv->push_back(1); // push(@pem,1);

			} else if (in.find("XC:i:") == string::npos && CIGAR.find("S") != string::npos) { // elsif($in!~/XC:i:/ and $CIAGR=~/S/){
				// $row="$list[0]";
				// for(my $j=1; $j<=8; $j++){
				// 	$row=$row."\t$list[$j]";
				// }
				// $row=$row."\tA\tH\tMD:Z:$CIAGR";
				row = list->at(0) + "\t" + list->at(1) + "\t" + list->at(2) + "\t" + list->at(3) + "\t" + list->at(4) + "\t" + list->at(5) + "\t" + list->at(6) + "\t" + list->at(7) + "\t" + list->at(8) + "\tA\tH\tMD:Z:" + CIGAR;
				rowv->push_back(row); // push(@row,$row);
				posv->push_back(pos); // push(@pos,$pos);
				pemv->push_back(0); // push(@pem,0);
			}

			#ifdef DEBUG_READ_INFO
				cout << "\n\n";
				cout << "in: " << in << "\n";
				cout << "pos: " << pos << "\n";
				cout << "mq: " << mq << "\n";
				cout << "CIGAR: " << CIGAR << "\n";
				cout << "ins_read: " << ins_read << "\n";
				cout << "bq: " << bq << "\n";
				cout << "alg_length: " << alg_length << "\n";
				cout << "ins_read: " << ins_read << "\n";
				cout << "ins_cut: " << ins_cut << "\n";
				cout << "large_pem: " << large_pem << "\n";
				cout << "row: " << row << "\n";
			#endif

			// while($CIAGR=~/(\d+)[MI]/g){
				// $alg_length+=$1;
			// }
			alg_length = count_CIGAR(CIGAR, true);

			#ifdef DEBUG_READ_INFO
				cout << "New alg_length: " << alg_length << "\n";
			#endif
			vector<string>().swap(*list);
			delete list;
			rc(pos, alg_length, mq, win_size_del); // rc($pos,$alg_length,$mq);
		}
	}
	
	pclose(sam_stream); // close(ST) ||die"$!\n";
		
	for (int j = 0; j < rowv->size(); j++) { // for(my $j=0; $j<=$#row; $j++){
		
		int last;
		int over = 0;

		if (j == 0) {
			last = rowv->size() - 1;
		} else {
			last = j - 1;
		}

		if (j != rowv->size() - 1) {
			over = posv->at(j+1);
		}
		
		if (pemv->at(j) == 1 || posv->at(j) - posv->at(last) < read_length || over - posv->at(j) < read_length) {
			PSAM << rowv->at(j) << "\n";
		}
	}
	PSAM.close(); // close(PSAM) or die "Can't close $chr_pem_sam";

	int start, end; // my ($start, $end);

	ofstream RC;
	RC.open(chr_rc, ios::out | ios::app); // open(RC,">$chr_rc") or die "Can't open $chr_rc";

	double mq_avg;
	double rc_avg;

	for (int i = 0; i < rcv->size(); i++) { // for(my $i=0; $i<=$#rc; $i++){
		mq_avg = 60; // $mq_avg=60;
		if (mqv->at(i) != -1) { // if(defined($mq[$i])){
			mq_avg = (double) mqv->at(i) / rcv->at(i); // $mq_avg=$mq[$i]/$rc[$i];
		}

		vector<int>* rd_window = new vector<int>(); // my @rd_window;
		start = i*win_size_del + 1; // $start=$i*$win_size_del+1;
		end = (i + 1) * win_size_del; // $end=($i+1)*$win_size_del;

		for (int j = start; j <= end && j < rd_bp->size(); j++) { // for(my $i=$start; $i<=$end; $i++){
			rd_window->push_back(rd_bp->at(j)); // push(@rd_window, $rd_bp[$i]);
		}

		rc_avg = (double) (midhalf(rd_window) / win_size_del); // $rc_avg=midhalf(@rd_window)/$win_size_del;

		char *buffer = new char[1024];

		sprintf(buffer, "%d\t%d\t%.2f\t%d\n", start, end, rc_avg, (int) (mq_avg + 0.5)); // printf RC ("$start\t$end\t%.2f\t%d\n",$rc_avg,int($mq_avg+0.5));
		RC << buffer;
		vector<int>().swap(*rd_window);
		delete rd_window;
	}

	RC.close(); // close(RC) or die "Can't close $chr_rc";
	

	return 0;
}
