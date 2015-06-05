<?php
// Modified by kerker!

/*
$tempf = file_get_contents("/input");
$lines = explode("\n","$tmepf");
$fi=fopen("input.php","w");
fputs($fi,"<?\n".'$seq1="'.$line[0].'";'."\n");
fputs($fi,'$seq2="'.$line[1].'";'."\n");
fputs($fi,'$matfile="../'.$line[2].'";'."\n");
fputs($fi,'$opp='.$line[3].'";'."\n");
fputs($fi,'$exp='.$line[4].'";'."\n");
fputs($fi,'$suboptimal='.$line[5].'";');
fclose($fi);

/*
 * function: semiglobal alignment
 * Input : ¤u§@¸ô®|
 * Output: result.php
 */

if(count($argv)!=2)
    {print "path??\n";exit();}
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
//print "$matfile";
$matrix=readmat($matfile);

$A=array();
$B=array();
$C=array();
//print "$seq1\n"; print "$seq2\n";
//for first "-" (row)
for($i=1;$i<=strlen($seq2);$i++){
    $A[0][$i]=-1000000; 
	$B[0][$i]=-1000000; 
	$C[0][$i]=-0;
	$score[0][$i]=0;
    $dir[0][$i]="C";
	$dirA[0][$i]="C";
	$dirB[0][$i]="C";
	$dirC[0][$i]="C";
}
//for first "|" (column)
for($j=1;$j<=strlen($seq1);$j++){
    $A[$j][0]=-1000000; 
	$B[$j][0]=0; 
	$C[$j][0]=-1000000;
	$score[$j][0]=0;
    $dir[$j][0]="B";
	$dirA[$j][0]="B";
	$dirB[$j][0]="B";
	$dirC[$j][0]="B";

}

$A[0][0]=0;
$B[0][0]=0;
$C[0][0]=0;
$score[0][0]=0;
$dir[0][0]="o";

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
if($total_len > 4000 && file_exists("$prefix/dirA_semi")) {@unlink("$prefix/dirA_semi");@unlink("$prefix/dirB_semi");@unlink("$prefix/dirC_semi");@unlink("$prefix/dir_semi");@unlink("$prefix/score_semi");}
// Table A B C
for($j=1;$j<=strlen($seq1);$j++){
    for($i=1;$i<=strlen($seq2);$i++){
        $x=array_key_exists($seq1[$j-1],$matrix)?$seq1[$j-1]:"*";
        $y=array_key_exists($seq2[$i-1],$matrix)?$seq2[$i-1]:"*";
	//	print $j."  ".$i."  ".$matrix[$x][$y]."\n";

        $temp=getmax(array($A[$j-1][$i-1],$B[$j-1][$i-1],$C[$j-1][$i-1]));
        $A[$j][$i]=$matrix[$x][$y]+$temp;
        $As=$A[$j-1][$i-1];
        $Bs=$B[$j-1][$i-1];
        $Cs=$C[$j-1][$i-1];
		$dirA_t=determingABC($As,$Bs,$Cs,$temp);
		if($total_len>4000){
			file_put_contents("$prefix/dirA_semi",$j."\t".$i."\t".$dirA_t."\n",FILE_APPEND);
		}
 		else
			$dirA[$j][$i] = $dirA_t;

		$temp=getmax(array($A[$j-1][$i]+$opp+$exp,$B[$j-1][$i]+$exp,$C[$j-1][$i]+$opp+$exp));
        $B[$j][$i]=$temp;
        $As=$A[$j-1][$i]+$opp+$exp;
        $Bs=$B[$j-1][$i]+$exp;
        $Cs=$C[$j-1][$i]+$opp+$exp;
		$dirB_t=determingABC($As,$Bs,$Cs,$temp);
		if($total_len>4000){
			file_put_contents("$prefix/dirB_semi",$j."\t".$i."\t".$dirB_t."\n",FILE_APPEND);
		}
 		else
			$dirB[$j][$i] = $dirB_t;

		$temp=getmax(array($A[$j][$i-1]+$opp+$exp,$B[$j][$i-1]+$opp+$exp,$C[$j][$i-1]+$exp));
        $C[$j][$i]=$temp;
        $As=$A[$j][$i-1]+$opp+$exp;
        $Bs=$B[$j][$i-1]+$opp+$exp;
        $Cs=$C[$j][$i-1]+$exp;
		$dirC_t=determingABC($As,$Bs,$Cs,$temp);
		if($total_len>4000){
			file_put_contents("$prefix/dirC_semi",$j."\t".$i."\t".$dirC_t."\n",FILE_APPEND);
		}
 		else
			$dirC[$j][$i] = $dirC_t;

		$score[$j][$i]=getmax(array($A[$j][$i],$B[$j][$i],$C[$j][$i]));
        if($C[$j][$i]==$score[$j][$i])  $dir_t="C";//.$j.":".$i.":".$dirC[$j][$i];
        if($B[$j][$i]==$score[$j][$i])  $dir_t="B";//.$j.":".$i.":".$dirB[$j][$i];
        if($score[$j][$i]==$A[$j][$i])  $dir_t="A";//.$j.":".$i.":".$dirA[$j][$i];
		if($total_len>4000){
			file_put_contents("$prefix/dir_semi",$j."\t".$i."\t".$dir_t."\n",FILE_APPEND);
		}
 		else{
			$dir[$j][$i] = $dir_t;
		}
	}
}
$A = NULL;
$B = NULL;
$C = NULL;
unset($A);
unset($B);
unset($C);
//print_r($score);
//print_r($dir);
/*
for($i=1;$i<=strlen($seq2);$i++){
    print "\n";
    for($j=1;$j<=strlen($seq1);$j++){
        print $A[$j][$i]."  ";
    }
}
print "\n";
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
        print $B[$j][$i]."  ";
    }
}
print "\n";for($i=1;$i<=strlen($seq2);$i++){
    print "\n";
    for($j=1;$j<=strlen($seq1);$j++){
        print $dirB[$j][$i]."  ";
    }
}
print "\n";

for($i=1;$i<=strlen($seq2);$i++){
    print "\n";
    for($j=1;$j<=strlen($seq1);$j++){
        print $C[$j][$i]."  ";
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


//print the score and direction results		
for($i=0;$i<=strlen($seq2);$i++){
    print "\n";
    for($j=0;$j<=strlen($seq1);$j++){
        print $score[$j][$i]."  ";
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


//sorting by the score
/*
if($total_len > 4000){
	$data = file_get_contents("score_semi");
	$lines = explode("\n",$data);
	foreach($lines as $row){
		if($row == "") continue;
		$each = explode("\t",$row);
		$score[$each[0]][$each[1]] = $each[2];
	}
}	
*/
$set=array();
for($i=1;$i<=strlen($seq2)-1;$i++){	$set["v$i"]=$score[strlen($seq1)][$i]; } //from last column (|)
for($i=1;$i<=strlen($seq1);$i++){	$set["h$i"]=$score[$i][strlen($seq2)]; } //from last row (-)
/*
if($total_len > 4000){
	$score = NULL;
	unset($score);
}
*/
//print_r($set);
arsort($set,SORT_NUMERIC);
if($total_len > 4000){
	$file = file_get_contents("$prefix/dir_semi");
	$lines = explode("\n",$file);
	foreach($lines as $row){
		if($row == "") continue;
		$each = explode("\t",$row);
		$dir[$each[0]][$each[1]] = $each[2];
	}
	$file = NULL;
	$lines = NULL;
	unset($file);
	unset($lines);
}

$Cset=array();
//print_r($set);
foreach($set as $k=>$s){
	if(count($Cset)>=$suboptimal && $s!=$tempscore) break;
	if($k[0]=="h"){ //for last row (-)
		$i=strlen($seq2); $j=substr($k,1);
	}
	else{ //for last column (|)
		$i=substr($k,1); $j=strlen($seq1);
	}
	//print "$j:$i\n";
	if(($k[0]=="h"&& substr($dir[$j][$i],0,1)=="B")||($k[0]=="v"&&substr($dir[$j][$i],0,1)=="C"))	continue;
    array_push($Cset,$k);
	$tempscore=$s;
}
//print_r($Cset);
$set = NULL;
unset($set);
if($total_len > 4000){
	$dir = NULL;
	unset($dir);
}

//print_r($set);
//print_r($dir);

$result = array();
foreach($Cset as $k){
	if($k[0]=="v") {
    	$maxI=substr($k,1);
	    $maxJ=strlen($seq1);
	}
	else {
    	$maxI=strlen($seq2);
   		$maxJ=substr($k,1);
	}
	$i=$maxI; //print "i=$i ";
	$j=$maxJ; //print "j=$j\n";
	if($k[0]=="v"){ //for last column (|)
		$res2=substr($seq2,$maxI);
		$res1=str_repeat("-",strlen($seq2)-$maxI);
	}
	else{ //for last row (-)
		$res1=substr($seq1,$maxJ);
		$res2=str_repeat("-",strlen($seq1)-$maxJ);
	}
	//$lasti=$maxI;
	//$lastj=$maxJ;
	//print "$seq1\n";
	//print "$seq2\n";
//	$D = NULL;
//	unset($D);
	if($total_len > 4000){
		$D = NULL;
		unset($D);
		$D[0][0]="o";
		for($a=1;$a<=strlen($seq2);$a++){
			$D[0][$a] = "C";
		}
		for($b=1;$b<=strlen($seq1);$b++){
			$D[$b][0] = "B";
		}
		$dfile = fopen("$prefix/dir_semi","r");
		while($lines = fgets($dfile,64)){
			if($lines == "") continue;
			$each = explode("\t",$lines);
			$D[$each[0]][$each[1]] = trim($each[2]);
		}
		fclose($dfile);
	}
	else
		$D = $dir;
//	print_r($D);
//	if(isset($D)) echo "dir set!\n";
//	$D = $dir;
//	echo "D=".$D[$j][$i]."\n";
//	$dir = NULL;
//	unset($dir);
	$flag = "";	
	while($j+$i!=0){
//		echo "D=".$D[$j][$i]."j=".$j."i=".$i."\n";
		//**A**
		if($D[$j][$i]=="A"){
			if($total_len > 4000 && $flag != "A"){
				$dirB = NULL;
				$dirC = NULL;
				unset($dirB);
				unset($dirC);

				for($a=1;$a<=strlen($seq2);$a++){
					$dirA[0][$a] = "C";
				}
				for($b=1;$b<=strlen($seq1);$b++){
					$dirA[$b][0] = "B";
				}
				$data = file_get_contents("$prefix/dirA_semi");
				$lines = explode("\n",$data);
				foreach($lines as $row){
					if($row == "") continue;
					$each = explode("\t",$row);
					$dirA[$each[0]][$each[1]] = $each[2];
				}
				
			}
			if($dirA[$j][$i]==5 || $dirA[$j][$i]==1 || $dirA[$j][$i]==2 || $dirA[$j][$i]==3 ){
				$D[$j-1][$i-1]="A";
			}
			if($dirA[$j][$i]==6 || $dirA[$j][$i]==4 && $j>1 && $i>1){
                $D[$j-1][$i-1]="B";
            }
			if($dirA[$j][$i]==7 && $j>1 && $i>1){
                $D[$j-1][$i-1]="C";
            }
			$res1=$seq1[$j-1].$res1; $j--;
			$res2=$seq2[$i-1].$res2; $i--;
			
			$flag = "A";
	//		echo "D=".$D[$j-1][$i-1]."\n";
		}
		//**B**
		if($D[$j][$i]=="B"){
			if($total_len > 4000 && $flag != "B"){
				$dirA = NULL;
				$dirC = NULL; 
				unset($dirA);
				unset($dirC);
				for($a=1;$a<=strlen($seq2);$a++){
					$dirB[0][$a] = "C";
				}
				for($b=1;$b<=strlen($seq1);$b++){
					$dirB[$b][0] = "B";
				}
				$data = file_get_contents("$prefix/dirB_semi");
				$lines = explode("\n",$data);
				foreach($lines as $row){
					if($row == "") continue;
					$each = explode("\t",$row);
					$dirB[$each[0]][$each[1]] = $each[2];
				}
			}
			if($dirB[$j][$i]==5 || $dirB[$j][$i]==1 || $dirB[$j][$i]==2 || $dirB[$j][$i]==3 ){
                $D[$j-1][$i]="A";
            }
			if($dirB[$j][$i]==6 || $dirB[$j][$i]==4 && $j>1){
				$D[$j-1][$i]="B";
			}
			if($dirB[$j][$i]==7 && $j>1){
                $D[$j-1][$i]="C";
            }
			$res1=$seq1[$j-1].$res1; $j--;
			$res2="-".$res2;

			$flag = "B";
        }
		//**C**
		if($D[$j][$i]=="C"){
			if($total_len > 4000 && $flag != "C"){
				$dirA = NULL;
				$dirB = NULL;
				unset($dirA);
				unset($dirB);
				for($a=1;$a<=strlen($seq2);$a++){
					$dirC[0][$a] = "C";
				}
				for($b=1;$b<=strlen($seq1);$b++){
					$dirC[$b][0] = "B";
				}
				$data = file_get_contents("$prefix/dirC_semi");
				$lines = explode("\n",$data);
				foreach($lines as $row){
					if($row == "") continue;
					$each = explode("\t",$row);
					$dirC[$each[0]][$each[1]] = $each[2];
				}
			}
            if($dirC[$j][$i]==5 || $dirC[$j][$i]==1 || $dirC[$j][$i]==2 || $dirC[$j][$i]==3 ){
                $D[$j][$i-1]="A";
            }
            if($dirC[$j][$i]==6 || $dirC[$j][$i]==4 && $i>1){
                $D[$j][$i-1]="B";
            }
            if($dirC[$j][$i]==7 && $i>1){
                $D[$j][$i-1]="C";
            }
            $res1="-".$res1;
			$res2=$seq2[$i-1].$res2; $i--;

			$flag = "C";
        }
	}
//print $score[$maxJ][$maxI]."\n";
//print "$k \t$res1\n";
//print "$k\t$res2\n";
/*
	if($total_len > 4000){
		$score[0][0]=0;
		for($a=1;$a<=strlen($seq2);$a++){
			$score[0][$a] = "0";
		}
		for($b=1;$b<=strlen($seq1);$b++){
			$score[$b][0] = "0";
		}
		$data = file_get_contents("score_semi");
		$lines = explode("\n",$data);
		foreach($lines as $row){
			if($row == "") continue;
			$each = explode("\t",$row);
			$score[$each[0]][$each[1]] = $each[2];
		}
	}	
*/
	$Nmat = match_number($res1,$res2);	
	array_push($result,array("seq1"=>$res1,"seq2"=>$res2,"start1"=>0,"start2"=>0,"end1"=>strlen($seq1)-1,"end2"=>strlen($seq2)-1,"Nmat"=>$Nmat,"score"=>$score[$maxJ][$maxI]));
/*
	if($tatal_len > 4000){
		$score = NULL;
		unset($score);
	}
*/
}

$file = fopen($prefix."result.php","w");
//print $prefix."result.php";
if($file){
	fputs($file,gen_resultlist($result));
	fclose($file);
	}
else
	print "\nError when writing results!\n";

?>
