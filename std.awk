#std.awk
#usage:cat go-basic.obo | grep -Po "(?<!_)id: GO:[0-9]{7}|is_a: GO:[0-9]{7}" | awk -f std.awk >data_extract

BEGIN{
FS=":"
id=0
}

{
 if($1=="id"){
	if(id!=0&&length(is)!=0){
		print id
		for(item in is) printf("%s\t",is[item])
        printf("\n")
	delete is
    }
	flag=1
	id=$3
 }
 if($1=="is_a") is[flag++]=$3
}

END{
	print id
	for(item in is) printf("%s\t",is[item])
    printf("\n")
}
