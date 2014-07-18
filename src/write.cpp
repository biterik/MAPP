#include "write.h"
using namespace std;
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Write::Write(MAPP* mapp):InitPtrs(mapp)
{
    last_write_step=-1;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Write::~Write()
{
    
}
/*--------------------------------------------
 init before a run
 --------------------------------------------*/
void Write::init()
{
    if(last_write_step!=step_no)
    {
        write_file(step_no);
        last_write_step=step_no;
    }
    write_step=step_no+write_step_tally;
}
/*--------------------------------------------
 write the file
 --------------------------------------------*/
void Write::write()
{
    if(write_step!=step_no)
        return;
    write_file(step_no);
    last_write_step=step_no;
    write_step=step_no+write_step_tally;
}
/*--------------------------------------------
 after a run
 --------------------------------------------*/
void Write::fin()
{
    if(last_write_step!=step_no)
        write_file(step_no);
    last_write_step=step_no;
}