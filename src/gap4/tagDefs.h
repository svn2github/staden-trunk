#ifndef _tagDefs_h
#define _tagDefs_h

#define TAG_UNCHANGED         (0)
#define TAG_POSITION_CHANGED  (1L<<1)
#define TAG_LENGTH_CHANGED    (1L<<2)
#define TAG_TYPE_CHANGED      (1L<<3)
#define TAG_COMMENT_CHANGED   (1L<<4)
#define TAG_INSERTED          (1L<<5)
#define TAG_COMMENT_IN_MEMORY (1L<<7)
#define TAG_NEXT_DELETED      (1L<<8)

#endif  /*_tagDefs_h*/
