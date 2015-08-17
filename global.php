<?php

/*
 * function: semiglobal alignment
 * Input : Working Path
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


/*--------------資料是include就已經將變數宣告完成----------*/
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

/*---------------------------Initalize first column and first row's score-------------------*/
$A[0][0]=0;
$B[0][0]=0;
$C[0][0]=0;
$score[0][0]=0;
$dir[0][0]="o";

echo "\n\n".strlen($seq2)."  ".strlen($seq1)."\n";

//for first "-" (row)
for($i=1;$i<=strlen($seq2);$i++){
    $A[0][$i]=-1000000; 
	$B[0][$i]=-1000000; 
	$C[0][$i]=$opp+$exp*$i;
	$score[0][$i]=$opp+$exp*$i;
    $dir[0][$i]="C";
    $dirA[0][$i]="C";
    $dirB[0][$i]="C";
    $dirC[0][$i]="C";
}
//for first "|" (column)
for($j=1;$j<=strlen($seq1);$j++){
    $A[$j][0]=-1000000;
	$B[$j][0]=$opp+$exp*$j; 
	$C[$j][0]=-1000000;
	$score[$j][0]=$opp+$exp*$j;
    $dir[$j][0]="B";
    $dirA[$j][0]="B";
    $dirB[$j][0]="B";
    $dirC[$j][0]="B";
//echo var_dump($A);
}

echo "Tracing.......\n";
echo var_dump($A);
/*-------------------------------------------------*/
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
echo $prefix."!!!!!!";
if($total_len > 4000 && file_exists("$prefix/dirA_glo")) 
{
@unlink("$prefix/dirA_glo");
@unlink("$prefix/dirB_glo");
@unlink("$prefix/dirC_glo");
@unlink("$prefix/dir_glo");
@unlink("$prefix/score_glo");}

for($j=1;$j<=strlen($seq1);$j++){
    for($i=1;$i<=strlen($seq2);$i++){
        $x=array_key_exists($seq1[$j-1],$matrix)?$seq1[$j-1]:"*";
        $y=array_key_exists($seq2[$i-1],$matrix)?$seq2[$i-1]:"*";
		$temp11 = 	$matrix[$x][$y];
		//table A B C
		//print "$i   $j   $temp11\n";/*
        print var_dump($temp11);
        $temp=getmax(array($A[$j-1][$i-1],$B[$j-1][$i-1],$C[$j-1][$i-1]));
        $A[$j][$i]=$matrix[$x][$y]+$temp;
        $As=$A[$j-1][$i-1];
        $Bs=$B[$j-1][$i-1];
        $Cs=$C[$j-1][$i-1];
        $dirA_t=determingABC($As,$Bs,$Cs,$temp);
        /*
		if($total_len > 4000)
			file_put_contents("$prefix/dirA_glo",$j."\t".$i."\t".$dirA_t."\n",FILE_APPEND);
        else
         */
			$dirA[$j][$i] = $dirA_t;

        $temp=getmax(array($A[$j-1][$i]+$opp+$exp,$B[$j-1][$i]+$exp,$C[$j-1][$i]+$opp+$exp));
        $B[$j][$i]=$temp;
        $As=$A[$j-1][$i]+$opp+$exp;
        $Bs=$B[$j-1][$i]+$exp;
        $Cs=$C[$j-1][$i]+$opp+$exp;
        $dirB_t=determingABC($As,$Bs,$Cs,$temp);
        /*
        if($total_len > 4000)
			file_put_contents("$prefix/dirB_glo",$j."\t".$i."\t".$dirB_t."\n",FILE_APPEND);
        else
         */
			$dirB[$j][$i] = $dirB_t;
		
        $temp=getmax(array($A[$j][$i-1]+$opp+$exp,$B[$j][$i-1]+$opp+$exp,$C[$j][$i-1]+$exp));
        $C[$j][$i]=$temp;
        $As=$A[$j][$i-1]+$opp+$exp;
        $Bs=$B[$j][$i-1]+$opp+$exp;
        $Cs=$C[$j][$i-1]+$exp;
        $dirC_t=determingABC($As,$Bs,$Cs,$temp);
        /*
		if($total_len > 4000)
			file_put_contents("$prefix/dirC_glo",$j."\t".$i."\t".$dirC_t."\n",FILE_APPEND);
        else
         */
			$dirC[$j][$i] = $dirC_t;

        $score[$j][$i]=getmax(array($A[$j][$i],$B[$j][$i],$C[$j][$i]));
        if($C[$j][$i]==$score[$j][$i])  $dir_t="C";//.$j.":".$i.":".$dirC[$j][$i];
        if($B[$j][$i]==$score[$j][$i])  $dir_t="B";//.$j.":".$i.":".$dirB[$j][$i];
        if($score[$j][$i]==$A[$j][$i])  $dir_t="A";//.$j.":".$i.":".$dirA[$j][$i];
        /*
		if($total_len > 4000){
			file_put_contents("$prefix/dir_glo",$j."\t".$i."\t".$dir_t."\n",FILE_APPEND);
		//	file_put_contents("score_glo",$j."\t".$i."\t".$score_t."\n",FILE_APPEND);
        }
         */
		else{
			$dir[$j][$i] = $dir_t;
		//	$score[$j][$i] = $score_t;
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
print "\n"b
//printle($j!=0 && $i!=0){
$flag = "";
while($j+$i!=0){
    //**A**
//	echo "D=".$D[$j][$i]."i=$i"."j=$j\n";
    if($D[$j][$i]=="A"){
        /*
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
			$data = file_get_contents("$prefix/dirA_glo");
			$lines = explode("\n",$data);
			foreach($lines as $row){
				if($row == "") continue;
				$each = explode("\t",$row);
				$dirA[$each[0]][$each[1]] = $each[2];
			}
        }	
     */
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
    }
    //**B**
    if($D[$j][$i]=="B"){
        /*
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
			$data = file_get_contents("$prefix/dirB_glo");
			$lines = explode("\n",$data);
			foreach($lines as $row){
				if($row == "") continue;
				$each = explode("\t",$row);
				$dirB[$each[0]][$each[1]] = $each[2];
			}
        }
     */
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
        /*
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
			$data = file_get_contents("$prefix/dirC_glo");
			$lines = explode("\n",$data);
			foreach($lines as $row){
				if($row == "") continue;
				$each = explode("\t",$row);
				$dirC[$each[0]][$each[1]] = $each[2];
			}
        }
     */
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

//p(int "$k \t$res1\n";
//print "$k\t$res2\n";
/*
if($total_len > 4000){
	$data = file_get_contents("score_glo");
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

$file = fopen($prefix."result.php","w");
//print $prefix."result.php";
if($file){
	fputs($file,gen_resultlist($result));
	fclose($file);
	}
else
	print "\nError when writing results!\n";


?>
