#!/usr/bin/env bash

middleman build
rsync -ahvc --delete build/* gouda@blopker.com:public/dev/semiconductor
