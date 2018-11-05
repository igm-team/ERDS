sub printtime{
	my $printout=$_[0];
	my @p2ths=qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	my @time=localtime();
	my $year=1900+$time[5];
	if(length($time[0])==1){
		$time[0]="0$time[0]";
	}
	if(length($time[1])==1){
		$time[1]="0$time[1]";
	}
	my $time="$time[2]:$time[1]:$time[0] $p2ths[$time[4]] $time[3] $year";
	print "$time\t$printout\n"; 
}

sub date{
	my @time=localtime();
	$time[5]+=1900;
	$time[4]++;
	if(length($time[4])==1){
		$time[4]="0$time[4]";
	}
	if(length($time[3])==1){
		$time[3]="0$time[3]";
	}
	return "$time[5]$time[4]$time[3]"; 
}

sub mean {
	my (@array) = @_;
	my $count=$#array;
	if($count==-1){
		return 0;
		exit;
	}
	my $s=0;
	for(my $i=0; $i<=$count; $i++){
		$s+=$array[$i];
	}
	return $s/($count+1);
}

sub med { 
	my (@array) = @_; 
	my $count = $#array; 

	my @sarray = sort { $a <=> $b } @array; 
	if ($count % 2) { 
		return ($sarray[$count/2] + $sarray[$count/2 + 1]) / 2;  
	}
	else { 
		return $sarray[int($count/2)];
	} 
}
sub std_dev {
	my (@array) = @_;
	return $array[0] if ($#array==0);
	my $s=0;
	my $sq_s=0;
	for(my $i=0; $i<=$#array; $i++){
		$s+=$array[$i];
	}
	$s/=($#array+1);
	for(my $i=0; $i<=$#array; $i++){
		$sq_s+=($array[$i]-$s)*($array[$i]-$s);
	}
	return sqrt($sq_s/($#array));
}

sub bin {
	my $n=$_[0];
	my @a=(0,0,0,0,0,0,0,0);
	my $i=0;
	while($n>=2**($i)){
		$a[$i]=($n%(2**($i+1)))/(2**($i));
		$n-=$a[$i]*(2**($i));
		$i++;
	}
	return @a;
}

sub cluster{
	my ($start,$end,$count_thh,$dist_thh,@data)=@_;
	my $length=$end-$start+1;
	my (@left,@right);
	for(my $i=0; $i<=int($#data-1)/2; $i++){
		$left[$i]=$data[$i];
		$right[$i]=$data[$i+int($#data-1)/2+1];
	}
	my $num_pairs=$#left+1;
	my (@center_l,@center_r);
	my $num_clust=int($num_pairs/$count_thh)+1;
	for(my $i=0; $i<=$num_clust-1; $i++){
		push @center_l, $left[int($i*$num_pairs/$num_clust)];
		push @center_r, $right[int($i*$num_pairs/$num_clust)];
	}
	my (@cluster_final_l,@cluster_final_r);
	my $iter=0;
	while($iter<10){
		$iter++;
		my (@cluster_l,@cluster_r);
		for(my $i=0; $i<=$num_pairs-1; $i++){
			my $closest=0;
			my $dist=10000;
			for(my $j=0; $j<=$num_clust-1; $j++){
				my $temp_dist=abs($left[$i]-$center_l[$j])+abs($right[$i]-$center_r[$j]);
				if ($temp_dist<$dist) {
					$dist=$temp_dist;
					$closest=$j;
				}
			}
			push (@cluster_l, [$left[$i], $closest]);
			push (@cluster_r, [$right[$i], $closest]);			
		}
		for(my $j=0; $j<=$num_clust-1; $j++){
			my @members_l=grep {$_->[1]==$j} @cluster_l;
			my @members_r=grep {$_->[1]==$j} @cluster_r;
			my ($sum_l,$sum_r);
			for(my $k=0; $k<=$#members_l; $k++){
				$sum_l+=$members_l[$k][0];
				$sum_r+=$members_r[$k][0];
			}
			if($#members_l>=0){
				my $new_center_l=$sum_l/($#members_l+1);
				my $new_center_r=$sum_r/($#members_l+1);
				$center_l[$j]=$new_center_l;
				$center_r[$j]=$new_center_r;
			}
		}
		@cluster_final_l=@cluster_l;
		@cluster_final_r=@cluster_r;
	}
	
	my $max_members=$count_thh;
	my $max_overlap=0;
	my $max_index=-1;
	for(my $j=0; $j<=$num_clust-1; $j++){
		my @members_l=grep {$_->[1]==$j} @cluster_final_l;
		my @members_r=grep {$_->[1]==$j} @cluster_final_r;
		my $cluster_total=$#members_l+1;
		next if($cluster_total<$count_thh);
		for(my $k=0; $k<=$cluster_total-1; $k++){
			if($members_l[$k][0]-$center_l[$j]>$dist_thh or $members_r[$k][0]-$center_r[$j]>$dist_thh){
				splice(@members_l, $k, 1);
				splice(@members_r, $k, 1);
				$k--;
				$cluster_total--;
			}
		}
		next if($#members_l+1<$count_thh);
		my ($sum_l,$sum_r);
		for(my $k=0; $k<=$#members_l; $k++){
			$sum_l+=$members_l[$k][0];
			$sum_r+=$members_r[$k][0];
		}
		my $new_center_l=$sum_l/($#members_l+1);
		my $new_center_r=$sum_r/($#members_l+1);
		$center_l[$j]=$new_center_l;
		$center_r[$j]=$new_center_r;
		my $overlap=min($center_r[$j],$end)-max($center_l[$j],$start);
		my $rate=min($overlap/$length,$overlap/($center_r[$j]-$center_l[$j]+1));
		#if($#members_l+1>$max_members or ($#members_l+1==$max_members and $rate>=$max_overlap)){
		if($rate*($#members_l+1)>=$max_overlap*$max_members){
			$max_index=$j;
			$max_members=$#members_l+1;
			$max_overlap=$rate;
		}
	}
	return 0 if($max_index==-1);
	return (int($center_l[$max_index]+0.5),int($center_r[$max_index]+0.5));
}

sub sortfile { 
	my $file_tobe_sort=$_[0];
	my $temp="$file_tobe_sort.tmp";
	my @list;
	my $in;
	my %chr;
	for(my $j=0; $j<=$#chr; $j++){
		$chr{$chr[$j]}=$j;
	}

	open(F, $file_tobe_sort) or die "cannot open $file_tobe_sort";
	my @in;
	my @sin;
	my %hash1;
	my %hash2;
	my %hash3;
	while($in=<F>){
		$in=~s/[\r\n]//g;
		next if(length($in)==0);
		@list=split /\s/, $in;
		$hash1{$in}=$chr{$list[0]};
		$hash2{$in}=$list[1];
		$hash3{$in}=$list[2];
		push(@in,$in);
	}
	close(F);	

	open(T, ">$temp") or die "cannot open $temp";
	@sin= sort { $hash1{$a} <=> $hash1{$b} || $hash2{$a} <=> $hash2{$b} || $hash3{$a} <=> $hash3{$b}} @in;
	print T join ("\n",@sin);
	print T "\n";
	close(T);
	my $diff=open(DIFF, "diff $file_tobe_sort $temp |");
	$in=<DIFF>;
	if(defined($in)){
		system("rm $file_tobe_sort");
		system("mv $temp $file_tobe_sort");		
	}
	else{
		system("rm $temp");
	}
}

sub mergefile { 
	my $file_tobe_merge=$_[0];
	my $temp="$file_tobe_merge.tmp";
	my @list;
	my ($in, $start, $end);
	my $overlap_rate_thh1=0.5;
	my $overlap_rate_thh2=0.9;

	open(F, $file_tobe_merge) or die "cannot open $file_tobe_merge";
	my (@rows,@start,@end);
	while($in=<F>){
		@list=split /\s/, $in;
		$start=$list[1];
		$end=$list[2];
		push(@start,$start);
		push(@end,$end);
		push(@rows,"$in");
	}
	close(F);
	
	my $num_rows=$#start;
	for(my $j=0; $j<=$num_rows; $j++){
		my $overlap=min($end[$j], $end[$j+1])-max($start[$j], $start[$j+1]);
		if($overlap>0){
			my $lengtha=$end[$j]-$start[$j]+1;
			my $lengthb=$end[$j+1]-$start[$j+1]+1;
			if($overlap/$lengtha>=$overlap_rate_thh1 and $overlap/$lengthb>=$overlap_rate_thh1){
				if($log_sum[$j]>$log_sum[$j+1]){
					splice (@rows, $j+1,1);
					splice (@start, $j+1,1);
					splice (@end, $j+1,1);
					$j--;
					$num_rows--;
				}
				elsif($log_sum[$j]<=$log_sum[$j+1]){
					splice (@rows, $j,1);
					splice (@start, $j,1);
					splice (@end, $j,1);
					$j--;
					$num_rows--;
				}
				else{
					if($lengtha<$lengthb){
						splice (@rows, $j+1,1);
						splice (@start, $j+1,1);
						splice (@end, $j+1,1);
						$j--;
						$num_rows--;
					}
					else{
						splice (@rows, $j,1);
						splice (@start, $j,1);
						splice (@end, $j,1);
						$j--;
						$num_rows--;
					}				
				}
			}
			elsif($overlap/$lengtha>=$overlap_rate_thh2){
				splice (@rows, $j,1);
				splice (@start, $j,1);
				splice (@end, $j,1);
				$j--;
				$num_rows--;
			}
			elsif($overlap/$lengthb>=$overlap_rate_thh2){
				splice (@rows, $j+1,1);
				splice (@start, $j+1,1);
				splice (@end, $j+1,1);
				$j--;
				$num_rows--;
			}
		}
	}
	
	open(T, ">$temp") or die "cannot open $temp";
	for(my $j=0; $j<=$num_rows; $j++){
		print T "$rows[$j]";
	}
	close(T);
	my $diff=open(DIFF, "diff $file_tobe_merge $temp |");
	$in=<DIFF>;
	if(defined($in)){
		system("rm $file_tobe_merge");
		system("mv $temp $file_tobe_merge");		
	}
	else{
		system("rm $temp");
	}
}

1
