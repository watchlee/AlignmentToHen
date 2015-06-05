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
 * ��X�}�C���̤j��
 * Input: �Ʀr�}�C
 * return: �}�C���̤j��
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
 * ����result.php�ɮ�
 * Input: result�}�C���
 * return: ��Jresult.php�ɮת��r��
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
 * Ū�JSubstitution Matrix
 * Input: SM���|
 * return: �Halphabet�����ޭȪ�2D Array
 */
function readmat($file){
	
	$fp=fopen($file,"r");
	if($fp==false) return false;
	//�p�G�}���ɮץ��ѴN�^��false
	$matrix=array();
	$col=array();
	fgets($fp); 		//���L�Ĥ@�����
	$line=fgets($fp);
	$element = strtok(trim($line)," ");
	while($element!==false){
		if($element=="") {$element = strtok(" ");continue;}
		array_push($col,trim($element));
		$element = strtok(" ");
	}
	while($line){
		$rec=array();
		$element = strtok(trim($line)," ");
    	while($element!==false){
        	if($element=="") {$element = strtok(" ");continue;}
        	array_push($rec,trim($element));
        	$element = strtok(" ");
    	}
		for($i=1;$i<count($rec);$i++)
			$matrix[$rec[0]][$col[$i-1]]=$rec[$i];		
		$line=fgets($fp);
	}
	fclose($fp);
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
