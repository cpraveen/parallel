# Developed by : Pritam Giri
# Date : 8.10.2016
# TIFR-CAM
t=1;
flag=1;
while [ "$flag" = "1" ];
do
	if ls $(printf "Solution-%0.4d-*" $t) >/dev/null 2>&1;
	then
		cat $(printf "Solution-%0.4d-*" $t) > $(printf "Solution-%0.4d.tec" $t)
		rm $(printf "Solution-%0.4d-*" $t)
		t=$((t + 1))
	else
		flag=0
	fi
done
