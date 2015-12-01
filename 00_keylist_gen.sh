#!/bin/bash

#Usage="15_keylist_gen.sh fieldname pwd "

# Output is
#fnum exptime JD FIL LONG LAT RAc DECc HAc ZDc BJDc
#1    2       3  4   5    6   7   8    9   10  11
# 


# This is a hacked version to be used only for G101_PR3

if (( $# != 2 )) ; then
cat > /dev/stderr <<EOF
Usage: $0 fieldname projid
EOF
exit 1
fi

fieldname=$1
project_id=$2

PROM_DB_USER=hatmaster
PROM_DB_HOST=hat
PROM_DB_PWD=kortefa11
MAINPATH=/datab

MYSQL="mysql -u $PROM_DB_USER -p$PROM_DB_PWD -h $PROM_DB_HOST --skip-column-names --batch "

# Output is
#fnum exptime JD FIL LONG LAT RAc DECc HAc ZDc BJDc
#1    2       3  4   5    6   7   8    9   10  11
# 

$MYSQL -s -N -e "select \
	I.IMstid, I.IMfnum, I.IMtexp, J.JD, FL.FLname, S.SIlongitude, \
        S.SIlatitude, 15*A.rac, A.decc, A.hac, A.zc, J.BJD, I.IMjd \
	from (HATMASTER.Images = I, HATMASTER.Filters = FL, \
             HATMASTER.Sites = S, HATMASTER.Stations = ST, \
             HATCALIB.AstromPointing = A) left join HATCALIB.JulianDates = J \
        on (I.IMfnum = J.fnum and I.IMstid = J.station_id) \
	where I.IMstatus = 0 and \
	(I.IMobject like 'G%${fieldname}' || I.IMobject = '${fieldname}') and \
	I.IMimtype = 'object' and \
	I.IMfilid = FL.FLid and \
	I.IMstid = ST.STid and ST.STsiid =  S.SIid and \
        I.IMprid = ${project_id} and \
        I.IMfnum = A.fnum and I.IMstid = A.station_id \
	group by I.IMstid, I.IMfnum, I.IMccd \
	order by I.IMstid, I.IMfnum, I.IMccd;" | \
gawk '{print $1"-"$2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13;}' | \
while read name texp JD FLnam SIlongitude SIlatitude ra dec ha zc BJD IMJD ; do
     echo $JD | grep NULL > /dev/null && (
         pars=$(hjd 24$IMJD $ra $dec)
         echo $name $texp 24$IMJD $FLnam $SIlongitude $SIlatitude $ra $dec $ha $zc ${pars[0]}
     ) || (
         echo $name $texp 24$IMJD $FLnam $SIlongitude $SIlatitude $ra $dec $ha $zc $BJD
     )
done 
