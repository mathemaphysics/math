#ifndef SYMTAB_H
#define SYMTAB_H

#define SYMBOL_MAX_NAME_LEN 256 ///< The maximum number of bytes in symbol_t::name
#define SYMTAB_DEFAULT_SIZE 9997 ///< The fixed size of the symbol table used for configuring yamd

typedef struct symbol
{
    /**
     * The name of the variable. This is the character
     * string which is hashed to give the location of
     * the start of the bucket in which the values lie.
     */
    char name[SYMBOL_MAX_NAME_LEN]; ///< The variable name

    /**
     * This is the size in bytes of symbol_t::name.
     */
    long size; ///< The number of bytes to allocate for data

    /**
     * Points to the beginning of the data stored as
     * the value of the variable in the symbol table.
     */
    void *data; ///< Pointer to the data
} symbol_t;

typedef struct symtab
{
    /**
     * The total number of symbols allocated.
     */
    long size;

    /**
     * The total number of symbols currently occupying
     * the symbol table.
     */
    long nsyms;

    /**
     * This is the pointer to the beginning of the symbol
     * table itself.
     */
    symbol_t *syms; /// The start of the list of symbols

    /**
     * Each entry in symtab_t::occs contains either a 0
     * or a 1, indicating that the space is not occupied
     * and occupied, respectively.
     */
    char *occs; /// Occupancy indicator for each entry in symbol_t::syms
} symtab_t;

int symbol_init( symbol_t *, char *, int );

int symbol_free( symbol_t * );

unsigned int symbol_hash( char * );

int symbol_set_data( symbol_t *, void * );

int symtab_init( symtab_t *, long );

int symtab_free( symtab_t * );

int symtab_add( symtab_t *, char *, long, void * );

int symtab_del( symtab_t *, char * );

symbol_t *symtab_lookup( symtab_t *, char * );

void symtab_print( symtab_t * );

#endif

// vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab:smarttab
