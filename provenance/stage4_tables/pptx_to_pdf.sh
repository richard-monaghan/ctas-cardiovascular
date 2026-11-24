#!/bin/bash
# usage: pptx_to_pdf.sh <in.pptx> <out.pdf>  -- export via Microsoft PowerPoint (AppleScript); the input is a staged copy only
IN="$1"; OUT="$2"
osascript <<OSA
tell application "Microsoft PowerPoint"
  open (POSIX file "$IN")
  delay 2
  set p to active presentation
  save p in (POSIX file "$OUT") as save as PDF
  delay 1
  close p saving no
end tell
OSA
