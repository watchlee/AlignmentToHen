<?php
/*
 * function: semiglobal alignment
 * Input : path, T 
 * Output: result.php
 */

if(count($argv)!=3)
    {print "path?? or T?\n";exit();}
$prefix = $argv[1];
if($prefix[strlen($prefix)-1]!="/")
    $prefix.="/";

if(@fopen($prefix."input.php","r")===false)
{print "input file??\n";exit();}
include_once($prefix."input.php");
include_once("util.php");

/*
 * Input:
 * $seq1 $seq2: sequence with structural alphabet
 * $matfile: Substitution matrix file
 * $opp: open gap penelty
 * $exp: extended gap penalty
 * $suboptimal: # of output solution
*/
$matrix=readmat($matfile);
$As=array();
$Bs=array();
$Cs=array();

//**** fill in the initial values ****
//for first "-" (row)
for($i=1;$i<=strlen($seq2);$i++){
	for($k=-2;$k<=$i;$k++){
		if($k>0){
			$Al[0][$i][$k] = -1000000;
	        $Bl[0][$i][$k] = -1000000;
			$Cl[0][$i][$k] = $opp+$exp*$k;
		}	
		else{
			$Al[0][$i][$k]=0;
			$Bl[0][$i][$k]=0;
			$Cl[0][$i][$k]=0;
		}
    	$dir[0][$i][$k]="C";
		$dirA[0][$i][$k]="C";
    	$dirB[0][$i][$k]="C";
    	$dirC[0][$i][$k]="C";
	}
}
//for first "|" (column)
for($j=1;$j<=strlen($seq1);$j++){
	for($k=-2;$k<=$j;$k++){
		if($k>0){
			$Al[$j][0][$k] = -1000000;
			$Bl[$j][0][$k] = $opp+$exp*$k;
			$Cl[$j][0][$k] = -1000000;
		}
		else{
			$Al[$j][0][$k]=0;
			$Bl[$j][0][$k]=0;
			$Cl[$j][0][$k]=0;
		}
    	$dir[$j][0][$k]="B";
		$dirA[$j][0][$k]="B";
    	$dirB[$j][0][$k]="B";
    	$dirC[$j][0][$k]="B";
	}
}
$Al[0][0][0]=0;$Al[0][0][-1] = -1000000;//$Al[0][0][-2]=0;
$Bl[0][0][0]=0;$Bl[0][0][-1]= 0;//$Bl[0][0][-2]=$opp+$exp;
$Cl[0][0][0]=0;$Cl[0][0][-1]= 0;//$Cl[0][0][-2]=$opp+$exp;
$dir[0][0]="o";

//**** the function for determing the directions according where the score they got ****
function determingABC($As,$Bs,$Cs,$temp){
    if($Cs==$temp)  $determing="7";
    if($Bs==$temp)  $determing="6";
    if($As==$temp)  $determing="5";
    if($Cs==$temp && $Bs==$temp) $determing="4";
    if($As==$temp && $Cs==$temp) $determing="3";
    if($As==$temp && $Bs==$temp) $determing="2";
    if($As==$temp && $Bs==$temp && $Cs==$temp) $determing="1";
    return $determing;
}
$total_len = strlen($seq1)+strlen($seq2);
if($total_len > 300 && file_exists("$prefix/dirA_nlocal")) {@unlink("$prefix/dirA_nlocal");@unlink("$prefix/dirB_nlocal");@unlink("$prefix/dirC_nlocal");@unlink("$prefix/dir_nlocal");@unlink("$prefix/score_nlocal");}

//**** fill in the talbe A,B adn C ****
for($j=1;$j<=strlen($seq1);$j++){
    for($i=1;$i<=strlen($seq2);$i++){
        $x=array_key_exists($seq1[$j-1],$matrix)?$seq1[$j-1]:"*";
        $y=array_key_exists($seq2[$i-1],$matrix)?$seq2[$i-1]:"*";

		//**** table A ****
		$l=$j+$i;
		$temp=array();
		for($k=1;$k<=$l;$k++){
			if($k-2<0){
				$Al[$j-1][$i-1][$k-2] = -1000000;
				$Bl[$j-1][$i-1][$k-2] = -1000000;
				$Cl[$j-1][$i-1][$k-2] = -1000000;
			}
			$tempAl = $Al[$j-1][$i-1][$k-2]+$matrix[$x][$y];
            $tempBl = $Bl[$j-1][$i-1][$k-2]+$matrix[$x][$y];
            $tempCl = $Cl[$j-1][$i-1][$k-2]+$matrix[$x][$y];
			$Al[$j][$i][$k] = getmax(array($tempAl,$tempBl,$tempCl));
		//	$dirA[$j][$i][$k]=determingABC($tempAl,$tempBl,$tempCl,$Al[$j][$i][$k]);
			$dirA_t = determingABC($tempAl,$tempBl,$tempCl,$Al[$j][$i][$k]);
		
	        if($total_len>300){
		        file_put_contents("$prefix/dirA_nlocal",$j."\t".$i."\t".$k."\t".$dirA_t."\n",FILE_APPEND);
		    }
		    else
				$dirA[$j][$i][$k] = $dirA_t;
		}
		//**** table B ****
		$temp=array();
		for($k=1;$k<=$l;$k++){
			if($k-1==0){
				$Al[$j-1][$i][$k-1]=0;
				$Bl[$j-1][$i][$k-1]=$opp;
				$Cl[$j-1][$i][$k-1]=0;
			}
			$tempAl = $Al[$j-1][$i][$k-1]+$opp+$exp;
			$tempBl = $Bl[$j-1][$i][$k-1]+$exp;
			$tempCl = $Cl[$j-1][$i][$k-1]+$opp+$exp;
			$Bl[$j][$i][$k]=getmax(array($tempAl,$tempBl,$tempCl));
		//	$dirB[$j][$i][$k]=determingABC($tempAl,$tempBl,$tempCl,$Bl[$j][$i][$k]);
			$dirB_t = determingABC($tempAl,$tempBl,$tempCl,$Bl[$j][$i][$k]);
		
    	    if($total_len>300){
		        file_put_contents("$prefix/dirB_nlocal",$j."\t".$i."\t".$k."\t".$dirB_t."\n",FILE_APPEND);
		    }
		    else
				$dirB[$j][$i][$k] = $dirB_t;
		}
		//**** table C ****
		$temp=array();
		for($k=1;$k<=$l;$k++){
			if($k-1==0){
				$Al[$j][$i-1][$k-1]=0;
				$Bl[$j][$i-1][$k-1]=0;
				$Cl[$j][$i-1][$k-1]=$opp;
			}
			$tempAl = $Al[$j][$i-1][$k-1]+$opp+$exp;
            $tempBl = $Bl[$j][$i-1][$k-1]+$opp+$exp;
            $tempCl = $Cl[$j][$i-1][$k-1]+$exp;
			$Cl[$j][$i][$k]=getmax(array($tempAl,$tempBl,$tempCl));
		//	$dirC[$j][$i][$k]=determingABC($tempAl,$tempBl,$tempCl,$Cl[$j][$i][$k]);
			$dirC_t=determingABC($tempAl,$tempBl,$tempCl,$Cl[$j][$i][$k]);
		
        	if($total_len>300){
	        	file_put_contents("$prefix/dirC_nlocal",$j."\t".$i."\t".$k."\t".$dirC_t."\n",FILE_APPEND);
		    }
		    else
				$dirC[$j][$i][$k] = $dirC_t;
		}
		//**** decide the direction and record it ****
		$Sd=array();	
		for($k=1;$k<=$l;$k++){
			//print "$j:$i:$k\n";
			array_push($Sd,$Al[$j][$i][$k],$Bl[$j][$i][$k],$Cl[$j][$i][$k]);
			$scoredir[$j][$i][$k]=getmax($Sd);
	        if($Cl[$j][$i][$k] == $scoredir[$j][$i][$k]) $dir_t = "C";// $dir[$j][$i][$k]="C";
			if($Bl[$j][$i][$k] == $scoredir[$j][$i][$k]) $dir_t = "B";// $dir[$j][$i][$k]="B";
        	if($scoredir[$j][$i][$k] == $Al[$j][$i][$k]) $dir_t = "A";// $dir[$j][$i][$k]="A";		
        	if($total_len>300){
	        	file_put_contents("$prefix/dir_nlocal",$j."\t".$i."\t".$k."\t".$dir_t."\n",FILE_APPEND);
		    }
		    else
				$dir[$j][$i][$k] = $dir_t;
		}		
	}
}

//print_r($dirA);
//print_r($dirB);
//print_r($dirC);
//print_r($dir);
/*
//**** normalized the score of each candidate and put them into the array $tempset ****
//**** $tempset["index j":"index i":"length":"from which table"] = normalized score ****
//**** [Version1]: only eliminate the candidates with length shorter than input "T" ****
for($j=1;$j<=strlen($seq1);$j++){
    for($i=1;$i<=strlen($seq2);$i++){
		foreach($Al[$j][$i] as $k =>$stemp){
			if($k<$argv[2]) continue;
			//if($k<=0) continue;
			$NS=$stemp/$k;
			$tempset["$j:$i:$k:A"]=round($NS,2);
		}
		foreach($Bl[$j][$i] as $k =>$stemp){
			if($k<$argv[2]) continue;
			//if($k<=0) continue;
            $NS=$stemp/$k;
            $tempset["$j:$i:$k:B"]=round($NS,2);
        }
		foreach($Cl[$j][$i] as $k =>$stemp){
			if($k<$argv[2]) continue;
			//if($k<=0) continue;
            $NS=$stemp/$k;
            $tempset["$j:$i:$k:C"]=round($NS,2);
        }
	}
}
*/
//**** normalized the score of each candidate and put them into the array $tempset ****
//**** $tempset["index j":"index i":"length":"from which table"] = normalized score ****
//**** [Version2]: eliminate the repeat(shadow) candidate and the candidates with length shorter than input "T" ****
for($j=1;$j<=strlen($seq1);$j++){
    for($i=1;$i<=strlen($seq2);$i++){
		$set=array();
        foreach($Al[$j][$i] as $k =>$stemp){
            if($k<$argv[2]) continue;
            //if($k<=0) continue;
            $NS=$stemp/$k;
            $set["$j:$i:$k:A"]=round($NS,2);
        }
        foreach($Bl[$j][$i] as $k =>$stemp){
            if($k<$argv[2]) continue;
            //if($k<=0) continue;
            $NS=$stemp/$k;
            $set["$j:$i:$k:B"]=round($NS,2);
        }
        foreach($Cl[$j][$i] as $k =>$stemp){
            if($k<$argv[2]) continue;
            //if($k<=0) continue;
            $NS=$stemp/$k;
            $set["$j:$i:$k:C"]=round($NS,2);
        }
        if(isset($set)){
            arsort($set,SORT_NUMERIC);
            //print_r($tempset);
            foreach($set as $info=>$stemp){
                $tempset[$info]=$stemp;
                //print "$info = $stemp\n";
                break;
            }
        }
    }
}
$Al = NULL;
$Bl = NULL;
$Cl = NULL;
unset($Al);
unset($Bl);
unset($Cl);


//**** sortin the scores from biggest one to samllest
//print_r($tempset);
arsort($tempset,SORT_NUMERIC);
//print_r($tempset);

//**** print some information ****
/*
for($i=1;$i<=strlen($seq2);$i++){
    print "\n";
    for($j=1;$j<=strlen($seq1);$j++){
        print $As[$j][$i]."  ";
    }
}
print "\n";
for($i=1;$i<=strlen($seq2);$i++){
    print "\n";
    for($j=1;$j<=strlen($seq1);$j++){
        print $Bs[$j][$i]."  ";
    }
}
print "\n";
for($i=1;$i<=strlen($seq2);$i++){
    print "\n";
    for($j=1;$j<=strlen($seq1);$j++){
        print $Cs[$j][$i]."  ";
    }
}
print "\n";
//print_r($dir);
for($i=1;$i<=strlen($seq2);$i++){
    print "\n";
    for($j=1;$j<=strlen($seq1);$j++){
        print $dirA[$j][$i]."  ";
    }
}
print "\n";
for($i=1;$i<=strlen($seq2);$i++){
    print "\n";
    for($j=1;$j<=strlen($seq1);$j++){
        print $dirB[$j][$i]."  ";
    }
}
print "\n";
for($i=1;$i<=strlen($seq2);$i++){
    print "\n";
    for($j=1;$j<=strlen($seq1);$j++){
        print $dirC[$j][$i]."  ";
    }
}
print "\n";
*/

//**** print the score and direction results, but can't work in normalized local ****
//**** because it is an array in the position ****
/*
for($i=0;$i<=strlen($seq2);$i++){
    print "\n";
    for($j=0;$j<=strlen($seq1);$j++){
        print ($scoredir[$j][$i]."  ";
    }
}
print "\n";


for($i=0;$i<=strlen($seq2);$i++){
    print "\n";
    for($j=0;$j<=strlen($seq1);$j++){
        print $dir[$j][$i]."  ";
    }
}
print "\n";
*/
if($total_len > 300){
    $file = file_get_contents("$prefix/dir_nlocal");
    $lines = explode("\n",$file);
    foreach($lines as $row){
        if($row == "") continue;
        $each = explode("\t",$row);
        $dir[$each[0]][$each[1]][$each[2]] = $each[3];
    }
    $file = NULL;
    $lines = NULL;
    unset($file);
    unset($lines);
}


//**** backtracking from the tail of candidate sequence ****
$countsum=0;
$result = array();
foreach($tempset as $k=>$s){
	if($countsum >$suboptimal) break;
	$temp=explode(":",$k);
	$maxI=$temp[1];
	$maxJ=$temp[0];
	$L=	$temp[2];//**** the start length ****
	$i=$maxI;// print "i=$i ";
	$j=$maxJ;// print "j=$j\n";
	$res1="";
	$res2="";
	$dir[$j][$i][$L]=$temp[3];//**** the start direction **** 
	if($s<=0) continue;
	$flag = "";
	//print "i=$i j=$j L=$L\n";
	while($L>0){
		//**** direction A****
        if($dir[$j][$i][$L]=="A"){
			if($total_len > 300 && $flag != "A"){
				$dirB = NULL;
				$dirC = NULL;
				unset($dirB);
				unset($dirC);
				for($k=1;$k<=$l;$k++){
					for($a=1;$a<=strlen($seq2);$a++){
						$dirA[0][$a][$k] = "C";
					}
					for($b=1;$b<=strlen($seq1);$b++){
						$dirA[$b][0][$k] = "B";
					}
				}	
				$data = file_get_contents("$prefix/dirA_nlocal");
				$lines = explode("\n",$data);
				foreach($lines as $row){
					if($row == "") continue;
					$each = explode("\t",$row);
					$dirA[$each[0]][$each[1]][$each[2]] = $each[3];
				}
			}

            if($dirA[$j][$i][$L]==5 || $dirA[$j][$i][$L]==1 || $dirA[$j][$i][$L]==2 || $dirA[$j][$i][$L]==3 ){
                $dir[$j-1][$i-1][$L-2]="A";
            }
            if($dirA[$j][$i][$L]==6 || $dirA[$j][$i][$L]==4 && $j>1 && $i>1){
                $dir[$j-1][$i-1][$L-2]="B";
            }
            if($dirA[$j][$i][$L]==7 && $j>1 && $i>1){
                $dir[$j-1][$i-1][$L-2]="C";
            }
            $res1=$seq1[$j-1].$res1; $j--; $L--; if($L<=0){ $res2="-".$res2; break; }
            $res2=$seq2[$i-1].$res2; $i--; $L--; if($L<=0) break;
        }
        //**** direction B ****
        if($dir[$j][$i][$L]=="B"){
			if($total_len > 300 && $flag != "B"){
				$dirA = NULL;
				$dirC = NULL;
				unset($dirA);
				unset($dirC);
				for($k=1;$k<=$l;$k++){
					for($a=1;$a<=strlen($seq2);$a++){
						$dirB[0][$a][$k] = "C";
					}
					for($b=1;$b<=strlen($seq1);$b++){
						$dirB[$b][0][$k] = "B";
					}
				}	
				$data = file_get_contents("$prefix/dirB_nlocal");
				$lines = explode("\n",$data);
				foreach($lines as $row){
					if($row == "") continue;
					$each = explode("\t",$row);
					$dirB[$each[0]][$each[1]][$each[2]] = $each[3];
				}
			}
            if($dirB[$j][$i][$L]==5 || $dirB[$j][$i][$L]==1 || $dirB[$j][$i][$L]==2 || $dirB[$j][$i][$L]==3 ){
                $dir[$j-1][$i][$L-1]="A";
            }
            if($dirB[$j][$i][$L]==6 || $dirB[$j][$i][$L]==4 && $j>1){
                $dir[$j-1][$i][$L-1]="B";
            }
            if($dirB[$j][$i][$L]==7 && $j>1){
                $dir[$j-1][$i][$L-1]="C";
            }
            $res1=$seq1[$j-1].$res1; $j--; $L--; 
            $res2="-".$res2; if($L<=0) break;
        }
        //**** direction C ****
        if($dir[$j][$i][$L]=="C"){
			if($total_len > 300 && $flag != "C"){
				$dirB = NULL;
				$dirA = NULL;
				unset($dirB);
				unset($dirA);
				for($k=1;$k<=$l;$k++){
					for($a=1;$a<=strlen($seq2);$a++){
						$dirC[0][$a][$k] = "C";
					}
					for($b=1;$b<=strlen($seq1);$b++){
						$dirC[$b][0][$k] = "B";
					}
				}	
				$data = file_get_contents("$prefix/dirC_nlocal");
				$lines = explode("\n",$data);
				foreach($lines as $row){
					if($row == "") continue;
					$each = explode("\t",$row);
					$dirC[$each[0]][$each[1]][$each[2]] = $each[3];
				}
			}
            if($dirC[$j][$i][$L]==5 || $dirC[$j][$i][$L]==1 || $dirC[$j][$i][$L]==2 || $dirC[$j][$i][$L]==3 ){
                $dir[$j][$i-1][$L-1]="A";
            }
            if($dirC[$j][$i][$L]==6 || $dirC[$j][$i][$L]==4 && $i>1){
                $dir[$j][$i-1][$L-1]="B";
            }
            if($dirC[$j][$i][$L]==7 && $i>1){
                $dir[$j][$i-1][$L-1]="C";
            }
            $res1="-".$res1;
            $res2=$seq2[$i-1].$res2; $i--; $L--; if($L<=0) break;
        }
		
    }
	//print "$k \t$res1\n";
	//print "$k\t$res2\n";
	//**** count the length of the alignment region of seq1 and seq2 ****
	//**** don't count the gap ****
	$tempres1=""; $tempres2="";
	$temp=explode("-",$res1);
	foreach($temp as $add){	$tempres1 = $tempres1.$add;}
	$temp=explode("-",$res2);
    foreach($temp as $add){ $tempres2 = $tempres2.$add;}
	//	$tempres1=strtr($res1,"-","");	
	//	print "$tempres1\n";
	//	$tempres2=strtr($res2,"-","");  
	//	print "$tempres2\n";
	$lastJ=$j+strlen($tempres1)-1;	//print "$lastJ\n";
	$lastI=$i+strlen($tempres2)-1;	//print "$lastI\n";
	//	print "score=".round($score[$maxJ][$maxI]/$forT,3)."\n";
	//	print "res1=$res1\nres2=$res2\n";
	$Nmat = match_number($res1,$res2);
	array_push($result,array("seq1"=>$res1,"seq2"=>$res2,"start1"=>$j,"start2"=>$i,"end1"=>$lastJ,"end2"=>$lastI,"Nmat"=>$Nmat,"score"=>$s));
	$countsum++;
}

//**** open file and put in the results ****
$file = fopen($prefix."result.php","w");
//print $prefix."result.php";
if($file){
	fputs($file,gen_resultlist($result));
	fclose($file);
	}
else print "\nError when writing results!\n";

?>
