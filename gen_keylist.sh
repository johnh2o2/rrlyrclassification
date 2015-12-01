#!/bin/bash

# Purpose:
# Using MySQL and the 'hjd' program, generates a keylist for the field with
# crucial parameters, such as :

# SERIAL TEXP   JD              FIL LONG     LAT      RA         DEC      HA       ZD       HJD
#6-111679 300.00 2453763.9513028 I -110.87847 31.68142 160.066767 37.52714 1.272550 16.72601 2453763.9561222


Usage="gen_keylist.sh fieldname"
if [ $# -lt 1 ]; then
	echo $Usage
	exit 1
fi

PROM_DB1=HATMASTER
PROM_DB2=HAT

PROM_DB_USER1=hatuser
PROM_DB_USER2=hatread

fieldname=$1
PROM_DB_PWD=kortefa11
PROM_DB_HOST=hat

	if [[ $PROM_DB_HOST == "hat" ]] ; then
	    PROM_DB_USER=$PROM_DB_USER1
            PROM_DB=$PROM_DB1
	    RADECSTRING="A.IAra * 15.0, A.IAdec"
            RADECSTRING2=""
	else
	    PROM_DB_USER=$PROM_DB_USER2
            PROM_DB=$PROM_DB2
            RADECSTRING="coalesce(A.IAra * 15.0, A.IAnomra * 15.0), coalesce(A.IAdec,A.IAnomdec)"
            RADECSTRING2="and A.IAnomra is not null and A.IAnomdec is not null"
	fi

	MYSQL="mysql -u $PROM_DB_USER -p$PROM_DB_PWD -h $PROM_DB_HOST --skip-column-names --batch $PROM_DB "

	$MYSQL -s -N -e "select \
		I.IMstid,I.IMfnum, I.IMtexp, I.IMjd + 2400000, \
		FL.FLname, \
		S.SIlongitude, S.SIlatitude, \
		$RADECSTRING, A.IAHA, A.IAz \
		from Images = I, Filters = FL, ImAstrom = A, Sites = S, Stations = ST \
		where I.IMstatus = 0 and \
		(I.IMobject like 'G%${fieldname}' || I.IMobject = '${fieldname}') and \
		I.IMstid = A.IAstid and \
		I.IMfnum = A.IAfnum and \
		I.IMfilid = FL.FLid and \
		I.IMstid = ST.STid and ST.STsiid = S.SIid and \
		ST.STactive = 1 $RADECSTRING2 \
		group by IMstid, IMfnum \
		order by I.IMjd;" | \
gawk '{
	printf "%s-%s ",$1,$2;
	for (i=3; i<=NF; i++) {printf "%s ",$i;}
	printf "\n";
}' > $tempfile

cat $tempfile | \
gawk '{print $1,$3,$5,$6,$7,$8}' | \
while read id jd long lat ra dec; do
	pars=$(${HATPIPE}/bin/hjd $jd $ra $dec)
	echo $id ${pars[0]}
done | \
grmatch --match-id --col-ref-id 1 --col-inp-id 1 -r $tempfile -i - -o - | \
gawk '{
	for (i=1; i<=NF-2; i++) {printf "%s ",$i;}
	printf "%s\n", $NF;
}'

rm $tempfile
exit
