#ifndef CFGREAD_H
#define CFGREAD_H

#include "symtab.h"

int cfgread_load_symbols_f( char *, symtab_t * );

int cfgread_load_symbols_s( FILE *, symtab_t * );

int cfgread_save_symbols_f( char *, symtab_t * );

int cfgread_save_symbols_s( FILE *, symtab_t * );

#endif
