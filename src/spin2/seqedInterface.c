#include "seqedInterface.h"
#include "tkSeqedUtils.h"
#include "tkSeqed.h"

/*
 * Set the cursor position from an X,Y position within the sheet widget
 */
int seqedSetCursor(tkSeqed *se, int x, int y) {

#ifdef DEBUG
    printf("seqedSetCursor x %d y %d \n", x, y);
#endif
    if (y < 0 || y >= se->displayHeight || x < 0 || x >= se->displayWidth)
	return 1;
/*
    setCursorPosSeq(xx,
		    xx->displayPos - DB_RelPos(xx,seqList[y]) + x + 1,
		    seqList[y]);
*/
    seqed_setCursorPos(se, se->displayPos + x);

    if (se->cursorPos < se->extent_left)
	seqed_setCursorPos(se, se->extent_left);
    else
	if (se->cursorPos > se->extent_right)
	    seqed_setCursorPos(se, se->extent_right);

    /* 
     * always position cursor on row of sequence, irrespective of the actual
     * mouse position
     */
    seqed_positionCursor(se, se->cursorSeq, se->cursorPos);

    return 0;
}


/*
 * Move cursor right
 */
int seqedCursorRight(tkSeqed *se)
{
/*
    if (xx->editorState == StateDown)
	return 1;
*/
    /* check boundary cases */

    if ((se->cursorPos >= 1) && (se->cursorPos < se->extent_right)) {

	/* set cursor pos and do any callbacks */
	seqed_setCursorPos(se, se->cursorPos + 1);

	/* display cursor */
	seqed_showCursor(se, se->cursorSeq, se->cursorPos);
	return 0;
    } else {
	bell();
	return 1;
    }

}


/*
 * Move cursor left
 */
int seqedCursorLeft(tkSeqed *se)
{
/*
    if (xx->editorState == StateDown)
	return 1;
	*/

    if ((se->cursorPos > 1) && (se->cursorPos <= se->extent_right)) {
	seqed_setCursorPos(se, se->cursorPos - 1);
	seqed_showCursor(se, se->cursorSeq, se->cursorPos);
	return 0;
    } else {
	bell();
	return 1;
    }

}


/*
 * Move cursor down,
 * cycle if necessary
 */
int seqedCursorDown(tkSeqed *se)
{
/*    
    int *seqList,seqCount, cseq, cpos;
    int posInContig;
    int i;

    if (xx->editorState == StateDown)
	return 1;

    posInContig = positionInContig(xx,xx->cursorSeq,xx->cursorPos);
    seqList = sequencesInRegion(xx,posInContig-1,2);
    seqCount = linesInRegion(xx,posInContig-1,2);
    for(i=0;
	i<seqCount && seqList[i]!=xx->cursorSeq;
	i++);
    
    cseq = xx->cursorSeq;
    cpos = xx->cursorPos;
    caretDown2(xx, i, seqCount, &cseq, &cpos, seqList, posInContig);
    if (xx->cursorSeq != cseq || xx->cursorPos != cpos)
	setCursorPosSeq(xx, cpos, cseq);

    showCursor(xx, xx->cursorSeq, xx->cursorPos);
*/
    return 0;
}


/*
 * Move cursor up,
 * cycle if necessary
 */
int seqedCursorUp(tkSeqed *se)
{
/*
    int *seqList,seqCount, cseq, cpos;
    int posInContig;
    int i;

    if (xx->editorState == StateDown)
	return 1;

    posInContig = positionInContig(xx,xx->cursorSeq,xx->cursorPos);
    seqList = sequencesInRegion(xx,posInContig-1,2);
    seqCount = linesInRegion(xx,posInContig-1,2);
    for(i=0;
	i<seqCount && seqList[i]!=xx->cursorSeq;
	i++);
    
    cseq = xx->cursorSeq;
    cpos = xx->cursorPos;
    caretUp2(xx, i, seqCount, &cseq, &cpos, seqList, posInContig);
    if (xx->cursorSeq != cseq || xx->cursorPos != cpos)
	setCursorPosSeq(xx, cpos, cseq);

    showCursor(xx,xx->cursorSeq, xx->cursorPos);
*/
    return 0;
}    
