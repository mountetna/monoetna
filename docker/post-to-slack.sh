#!/bin/sh

set -e

NAME="$1"
CHANNEL="$2"
MSG="$3"
MSG=${MSG//\\/\\\\} # \ 
MSG=${MSG//\//\\\/} # / 
MSG=${MSG//\'/\\\'} # ' (not strictly needed ?)
MSG=${MSG//\"/\\\"} # " 
MSG=${MSG//   /\\t} # \t (tab)
MSG=${MSG//
/\\\n} # \n (newline)
MSG=${MSG//^M/\\\r} # \r (carriage return)
MSG=${MSG//^L/\\\f} # \f (form feed)
MSG=${MSG//^H/\\\b} # \b (backspace)

curl -X POST --data-urlencode "payload={\"channel\": \"#${CHANNEL}\", \"username\": \"$NAME\", \"text\": \"${MSG}\", \"icon_emoji\": \":ghost:\"}" "$SLACK_WEBHOOK_URL"

