#!/bin/sh
set -eu

DEST_DIR=${1:?destination directory is required}
TAX_URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
FETCHED_DATE="$(date -u +%F)"

mkdir -p "$DEST_DIR"

LAST_MODIFIED_RAW="$({ wget --server-response --spider -O /dev/null "$TAX_URL" 2>&1 || true; } | awk 'tolower($1) == "last-modified:" { sub(/^[^:]+: /, "", $0); sub(/\r$/, "", $0); print; exit }')"

if [ -n "$LAST_MODIFIED_RAW" ]; then
    TAX_DATE="$(date -u -d "$LAST_MODIFIED_RAW" +%F 2>/dev/null || printf '%s' "$FETCHED_DATE")"
else
    TAX_DATE="$FETCHED_DATE"
fi

wget -q "$TAX_URL" -O /tmp/taxdump.tar.gz

tar -xzf /tmp/taxdump.tar.gz -C "$DEST_DIR"
rm /tmp/taxdump.tar.gz

printf '%s\n' "$TAX_DATE" > "$DEST_DIR/VERSION"
printf 'source=%s\n' "$TAX_URL" >> "$DEST_DIR/VERSION"
printf 'last_modified=%s\n' "${LAST_MODIFIED_RAW:-unknown}" >> "$DEST_DIR/VERSION"
printf 'fetched=%s\n' "$FETCHED_DATE" >> "$DEST_DIR/VERSION"
