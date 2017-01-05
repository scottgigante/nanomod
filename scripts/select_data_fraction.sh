FRAC=$1
TRAIN_IN=$2
VAL_IN=$3

TRAIN_OUT="$TRAIN_IN.small"
VAL_OUT="$VAL_IN.small"

# count lines
TRAIN_N=$(cat $TRAIN_IN | wc -l)
VAL_N=$(cat $VAL_IN | wc -l)

# calculate fraction
TRAIN_FRAC=$(echo "($TRAIN_N * $FRAC)/1" | bc)
VAL_FRAC=$(echo "($VAL_N * $FRAC)/1" | bc)

# headers
head -n 1 $TRAIN_IN > $TRAIN_OUT
head -n 1 $VAL_IN > $VAL_OUT

# random selection
tail -n +2 $TRAIN_IN | shuf -n $TRAIN_FRAC - >> $TRAIN_OUT
tail -n +2 $VAL_IN | shuf -n $VAL_FRAC - >> $VAL_OUT
