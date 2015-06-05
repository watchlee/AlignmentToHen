<?php
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
//	print "As: $As, Bs: $Bs, Cs: $Cs temp: $temp\n";
	if($Cs==$temp)  $determing="7";
    if($Bs==$temp)  $determing="6";
    if($As==$temp)  $determing="5";
    if($Cs==$temp && $Bs==$temp) $determing="4";
    if($As==$temp && $Cs==$temp) $determing="3";
    if($As==$temp && $Bs==$temp) $determing="2";
    if($As==$temp && $Bs==$temp && $Cs==$temp) $determing="1";
//	print "determaing: $determing\n";
	return $determing;
}


for($j=1;$j<=strlen($seq1);$j++){
    for($i=1;$i<=strlen($seq2);$i++){
        $x=array_key_exists($seq1[$j-1],$matrix)?$seq1[$j-1]:"*";
        $y=array_key_exists($seq2[$i-1],$matrix)?$seq2[$i-1]:"*";
		$temp11 = $matrix[$x][$y];
		//table A B C
		//print "$i   $j   $temp11\n";/*
        $temp=getmax_local(array($A[$j-1][$i-1],$B[$j-1][$i-1],$C[$j-1][$i-1]));
        $A[$j][$i]=$matrix[$x][$y]+$temp;
        $As=$A[$j-1][$i-1];
        $Bs=$B[$j-1][$i-1];
        $Cs=$C[$j-1][$i-1];
        $dirA[$j][$i]=determingABC($As,$Bs,$Cs,$temp);

        $temp=getmax_local(array($A[$j-1][$i]+$opp+$exp,$B[$j-1][$i]+$exp,$C[$j-1][$i]+$opp+$exp));
        $B[$j][$i]=$temp;
        $As=$A[$j-1][$i]+$opp+$exp;
        $Bs=$B[$j-1][$i]+$exp;
        $Cs=$C[$j-1][$i]+$opp+$exp;
        $dirB[$j][$i]=determingABC($As,$Bs,$Cs,$temp);

        $temp=getmax_local(array($A[$j][$i-1]+$opp+$exp,$B[$j][$i-1]+$opp+$exp,$C[$j][$i-1]+$exp));
        $C[$j][$i]=$temp;
        $As=$A[$j][$i-1]+$opp+$exp;
        $Bs=$B[$j][$i-1]+$opp+$exp;
        $Cs=$C[$j][$i-1]+$exp;
        $dirC[$j][$i]=determingABC($As,$Bs,$Cs,$temp);

        $score[$j][$i]=getmax(array($A[$j][$i],$B[$j][$i],$C[$j][$i]));

        if($C[$j][$i]==$score[$j][$i])  $dir[$j][$i]="C";
        if($B[$j][$i]==$score[$j][$i])  $dir[$j][$i]="B";
        if($score[$j][$i]==$A[$j][$i])  $dir[$j][$i]="A";

		
	}
}

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
print "\n";
*/

//sorting by the score
$set=array();
for($i=1;$i<=strlen($seq2);$i++){
	for($j=1;$j<=strlen($seq1);$j++){
		$set["$j:$i"]=$score[$j][$i]; 
	}
} 
//print_r($set);
arsort($set,SORT_NUMERIC);
$Cset=array();
foreach($set as $k=>$s){
	if(count($Cset)>=$suboptimal) break;
    array_push($Cset,$k);
}
//print_r($Cset);

$result = array();
foreach($Cset as $k){
	$temp=explode(":",$k);
	$maxI=$temp[1];
	$maxJ=$temp[0];
	
	$i=$maxI; //print "i=$i ";
	$j=$maxJ; //print "j=$j\n";
	$res1="";
	$res2="";
	$D=$dir;
//	while($j!=0 && $i!=0 || $score[$j][$i]!=0){
	while($score[$j][$i]!=0){
        if($D[$j][$i]=="A"){
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
        }
        //**B**
        elseif($D[$j][$i]=="B"){
            if($dirB[$j][$i]==5 || $dirB[$j][$i]==1 || $dirB[$j][$i]==2 || $dirB[$j][$i]==3 ){
                $D[$j-1][$i]="A";
            }
            if($dirB[$j][$i]==6 || $dirB[$j][$i]==4 && $j>1){
                $D[$j-1][$i]="B";
            }
            if($dirB[$j][$i]==7 && $j>1){
                $D[$j-1][$i]="C";
            }
			//if($j=0) break;
            $res1=$seq1[$j-1].$res1; $j--;
            $res2="-".$res2;
        }
        //**C**
        elseif($D[$j][$i]=="C"){
            if($dirC[$j][$i]==5 || $dirC[$j][$i]==1 || $dirC[$j][$i]==2 || $dirC[$j][$i]==3 ){
                $D[$j][$i-1]="A";
            }
            if($dirC[$j][$i]==6 || $dirC[$j][$i]==4 && $i>1){
                $D[$j][$i-1]="B";
            }
            if($dirC[$j][$i]==7 && $i>1){
                $D[$j][$i-1]="C";
            }
			//if($i=0) break;
            $res1="-".$res1;
            $res2=$seq2[$i-1].$res2; $i--;
        }
    }
	if($res1[0]=="-"){	$res1=substr($res1,1,strlen($res1)-1); $res2=substr($res2,1,strlen($res2)-1); $i++; }
	if($res2[0]=="-"){	$res1=substr($res1,1,strlen($res1)-1); $res2=substr($res2,1,strlen($res2)-1); $j++;	}
//print "$k \t$res1\n";
//print "$k\t$res2\n";
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
	$Nmat = match_number($res1,$res2);
	array_push($result,array("seq1"=>$res1,"seq2"=>$res2,"start1"=>$j,"start2"=>$i,"end1"=>$lastJ,"end2"=>$lastI,"Nmat"=>$Nmat,"score"=>$score[$maxJ][$maxI]));
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
