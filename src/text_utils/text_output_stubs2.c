/*
 * This single routine is in a separate file simply because gap/xgap have
 * their own UpdateTextOutput (dialogues.c), but still need to link
 * the verror() function. Splitting the two resolves link conflicts.
 */
void UpdateTextOutput(void) {
}
