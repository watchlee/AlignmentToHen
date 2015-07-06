<?
function gapsa2nt($sa,$nt,$last){
    $flag=1;$r=$sa;
    for($i=0;$i<strlen($sa);$i++){
        if($sa[$i]=="-") continue;
        $r[$i]=$nt[$last-1+$flag++];
    }
    return $r;
}


/* function getmax
 * 找出陣列中最大值
 * Input: 數字陣列
 * return: 陣列中最大值
 */
function getmax($arr){

	foreach ($arr as $data){
		if(!isset($max)||$data>$max)
			$max=$data;
	}
	return isset($max)?$max:false;
}

function getmax_local($arr){
	foreach ($arr as $data){
		if(!isset($max)||$data>$max)
			$max=$data;
	}	
	if($max<0)
		$max = 0;
	return isset($max)?$max:false;
}
/* function gen_resultlist
 * 產生result.php檔案
 * Input: result陣列資料
 * return: 輸入result.php檔案的字串
 */

function gen_resultlist($result){
	$r='<?'."\n".'$result_list=array();';
	foreach($result as $i =>$item){
		foreach($item as $content =>$value){
			$r.="\n".'$result_list['.$i.']["'.$content.'"] ="'.$value.'";';

		}
	}
	$r.="\n?>";
	return $r;
}
/* function readmat
 * 讀入Substitution Matrix
 * Input: SM路徑
 * return: 以alphabet為索引值的2D Array
 */
function readmat($file){
	
	$fp=fopen($file,"r");
	if($fp==false) return false;
	//如果開啟檔案失敗就回傳false
	$matrix=array();
	$col=array();
    fgets($fp); 		//跳過第一行註解
    $line=fgets($fp);
    /*using trim function to delete meanlessing character, which was contained in input_string*/
    $element = strtok(trim($line)," ");
    /*Keeping take empty character until EOF*/
    /*Main goal is to put elements into array, and if there have any empty character or something character, which would cause array can't normal run, it should be token away!
     * get any elements*/

    /*This code is used to store first row line alphabet, i.e A B C E D F G H I J K L M N O P Q R S T X V W * */

	while($element!==false){
        if($element=="") {$element = strtok(" ");continue;}
		array_push($col,trim($element));
        $element = strtok(" ");
    }

    for($j = 0;$j < count($col);$j++)
        echo $col[$j];
    /*And this second code uses to first store alphabet, and then put argument into data array*/
    $count= 0;
    $loop = 0;
    while($line){
        $count = $count+1;
		$rec=array();
        $element = strtok(trim($line)," ");
        /*This code store column alphabet, cause it is a symmetric matrix*/
//        echo $element;
    	while($element!==false){
        	if($element=="") {$element = strtok(" ");continue;}
        	array_push($rec,trim($element));
            $element = strtok(" ");
        }
        /*This code store argument, which gets from input file*/
        for($i=1;$i<count($rec);$i++)
        {
            $loop++;
            /*Using row and column array to store data index value, etc: matrix[A][A] is store B,
             *but for C languge, we can't do such thing!
             *
             * */
            $matrix[$rec[0]][$col[$i-1]]=$rec[$i];	
     //       echo $rec[0]."and".$col[$i-1]."=".$rec[$i]." ";
//            echo " ".$rec[$i];
//            echo "array have count = ".count($rec);

        }
        echo $rec[0]."\n";
    //    echo "\n";
		$line=fgets($fp);
    }

    echo "\n";
//    echo $count." ".$loop;
    fclose($fp);
    //print $matrix[3][1];
	return $matrix;
}
function match_number($res1,$res2){
	$Nmat = 0;
	for($i=0;$i<strlen($res1);$i++){
		if($res1[$i] != "-" && $res2[$i] != "-"){
			$Nmat++;
		}
	}
	return $Nmat;
}
?>
