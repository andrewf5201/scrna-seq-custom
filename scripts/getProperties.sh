#!/bin/sh

function getProperty {
   PROP_KEY=$1
   PROPERTY_FILE=$2
   PROP_VALUE=`cat $PROPERTY_FILE | grep "^[^#]" | grep "$PROP_KEY" | cut -d'=' -f2`
   echo $PROP_VALUE
}
