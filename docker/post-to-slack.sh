#!/bin/bash

NAME="$1"
CHANNEL="$2"
MSG="$3"
MSG=${MSG//\\/\\\\}
MSG=${MSG//\//\\\/}
MSG=${MSG//\'/\\\'}
MSG=${MSG//\"/\\\"} 
MSG=${MSG//   /\\t}
MSG=${MSG//
/\\\n}
MSG=${MSG//^M/\\\r}
MSG=${MSG//^L/\\\f}
MSG=${MSG//^H/\\\b}

if ! [[ -z "$SLACK_WEBHOOK_URL" ]]; then
  curl -X POST --data-urlencode "payload={\"channel\": \"#${CHANNEL}\", \"username\": \"$NAME\", \"text\": \"${MSG}\", \"icon_emoji\": \":ghost:\"}" "$SLACK_WEBHOOK_URL"
fi

