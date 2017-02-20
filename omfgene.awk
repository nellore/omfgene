# Count discordant paired-end alignments & mates > 500k bases away

# Absolute value
function abs(v) {
	return v < 0 ? -v : v
}

# Round to nearest 100
function rnd(v) {
	r = v / 100
	i = int(r)
	d = r - i
	return d >= 0.5 ? (i + 1) * 100 : i * 100
}

# Parse SAM
# If primary alignment AND mate mapped AND (mate is on different chr OR mate is more than 500k away)
# print alignments rounded to nearest 100
# Format: mate 1 chrom TAB mate 1 rounded coordinate TAB mate 2 chrom TAB mate 2 rounded coordinate TAB 1 if secondary else 0
$7 != "*" && ($7 != "=" && $7 != $3 || (abs($4-$8) > 500000)) {
	print $3 "\t" rnd($4) "\t" ($7 == "=" ? $3 : $7) "\t" rnd($8) "\t" (int($2 / 256) % 2)
}
