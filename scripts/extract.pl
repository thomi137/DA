#!/usr/bin/perl 

use strict;
use warnings;

while(<>) {
     if(/^(((?:\d+?)\.(?:\d+)(?:e(?:-|)\d\d|)|\d+?(?:e(?:-|)\d\d|))\t((?:\d+?)\.(?:\d+)(?:e(?:-|)\d\d|)|\d+?(?:e(?:-|)\d\d|))\n)$/osm || /^(&\t&\n)$/osm || /^(&\n)$/osm){
        print (($1=~/&\n/osm) ? "&\t&\n" : $1);
    }
}



