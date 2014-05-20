
/* ============================================================================
 * p2/cell/obj.hh
 * ========================================================================= */

#ifndef objINCLUDED
#define objINCLUDED

#include "cell.hh"
#include <string>

/*
 * Return a cell for a given object file.
 * name -> the name of the object file to read the cell from
 * <- the cell corresponding to object file _name_
 */
Cell *objReadCell(const char *name);

/*
 * Return a cell for a given object content.
 * obj_file_contents -> contents of a .obj file read into a std::string
 * <- the cell corresponding to obj_file_contents
 */
Cell *objReadCellFromString(std::string obj_file_contents);

/*
 * Write a given cell to a given object file.
 * cell -> the cell to write
 * name  -> the name of the object file to write the cell to
 */
void objWriteCell(Cell *cell, const char *name);

/*
 * Return a copy of a given cell.
 * cell -> the cell to copy
 * <- a cell semantically equivalent to _cell_
 */
Cell *objCopyCell(Cell *cell);

/*
 * Return a copy of a given cell.
 * cell -> the cell to clone
 * <- a cell identical to _cell_
 */
Cell *objCloneCell(Cell *cell);

#endif /* #ifndef objINCLUDED */

