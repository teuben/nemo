// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2006                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================ 

#include <qmessagebox.h>
#include <iostream>
#include <qcheckbox.h>
#include <qcolordialog.h>
#include <qmessagebox.h>
#include <qfiledialog.h>
#include <qfileinfo.h>
#include <qtable.h>
#include <assert.h>

#define LOCAL_DEBUG 0
#include "print_debug.h"
const unsigned int MAX_OBJECTS=30;

// ============================================================================
void ParticlesSelectForm::init()
{
    table_range->setNumRows(MAX_OBJECTS);
    table_range->setColumnWidth( 0, 50 );
    table_range->setColumnWidth( 1, 80 );
    table_range->setColumnWidth( 2, 80 );
    table_range->setColumnWidth( 3, 40 );
    table_range->setColumnWidth( 4, 50 );

    table_list->setNumRows(MAX_OBJECTS);
    table_list->setColumnWidth( 0, 50 );
    table_list->setColumnWidth( 1, 160);
    table_list->setColumnWidth( 2, 40);
    table_list->setColumnWidth( 3, 50);
    
  // Connexion with the TABLE_RANGE  object
  connect( table_range, SIGNAL( clicked(int,int,int,const QPoint&) ),
	   this, SLOT( setColor(int,int) ) );
  connect( table_range, SIGNAL( currentChanged(int,int) ),
	   this, SLOT( currentChangedRange(int,int) ) );
  connect( table_range, SIGNAL( valueChanged(int,int) ),
	   this, SLOT( valueChanged(int,int) ) );
  
  // Connexion with the TABLE_LIST  object
  
  // file name dialog box if double clicked
  connect( table_list, SIGNAL( doubleClicked(int,int,int,const QPoint&) ),
	   this, SLOT( changeListFile(int,int) ) );
  // edit file name if clicked
  //connect( table_list, SIGNAL( clicked(int,int,int,const QPoint&) ),
//	   table_list, SLOT( editCell(int,int) ) );  
  
  connect( table_list, SIGNAL( currentChanged(int,int) ),
	   this, SLOT( currentChangedList(int,int) ) );
  // buttons
  connect( colorPushButton, SIGNAL( clicked() ), this, SLOT( setColor() ) );
  connect( okPushButton, SIGNAL( clicked() ), this, SLOT( accept()) );
  connect( applyPushButton, SIGNAL( clicked() ), this, SLOT( apply()) );
  connect( cancelPushButton, SIGNAL( clicked() ), this, SLOT( reject() ) );
   
}
// ============================================================================
void ParticlesSelectForm::destroy()
{
}

// ============================================================================
// SetParticlesRangeForm::accept()                                             
// accept [OK button]                                                          
void ParticlesSelectForm::accept()
{
  apply();
  QDialog::accept();
}

// ============================================================================
// SetParticlesRangeForm::apply()                                              
// collect all the data strored in the  tabls_range, fill particles range vector, and 
// send data to the glbox.                                                     
void ParticlesSelectForm::apply()
{
  my_psv->clear(); // Clear the vector
  VirtualParticlesSelect::nb_select=0;

  // proceed on range
  applyRange();
  // proceed on list
  applyList();
  emit applyData(&nbody,pos,my_psv);  // send signal to GLBox to update the data
}
// ============================================================================
// ParticlesSelectForm::applyRange()
// collect information from 'table_range' and fill my_psv accordingly
void ParticlesSelectForm::applyRange()
{
    for (unsigned int i=0; i < MAX_OBJECTS; i++) {
       bool ok1,ok2,commit=FALSE;
       int first = table_range->text( i, 1 ).toInt( &ok1 );
       int last  = table_range->text( i, 2 ).toInt( &ok2 );

       int step;
       if ( ok1 && ok2 ) {
          commit = TRUE;
          step= table_range->text( i, 3 ).toInt( &ok2 );
	  if (!ok2) step=1;    // force step to 1
       }	
    if (commit) {
      VirtualParticlesSelect * pr = new ParticlesRange();
      pr->first_part = first;
      pr->last_part  = last;
      pr->npart      = last+1-first;
      pr->step_part  = step;
      pr->defaultIndexTab();
      pr->col        = QColor( table_range->text( i, 4 ));
      pr->is_visible = ((QCheckBox*)table_range->cellWidget( i,0))->isChecked();
      ParticlesSelect * ps = new ParticlesSelect();
      ps->vps = pr;
      //pr->printRange();
      my_psv->push_back(*ps);
      PRINT_D cerr << "In SetParticlesRangeForm::apply(), my_psv->size() = "
                   <<  my_psv->size() << "\n";
      delete ps;
    }
  }
}
// ============================================================================
// ParticlesSelectForm::applyList()
// collect information from 'table_list' and fill my_psv accordingly
 
void ParticlesSelectForm::applyList()
{
    for (unsigned int i=0; i < MAX_OBJECTS; i++) {
	QString filename = table_list->text( i,1); 
	if (filename!="") {	// check if a filename exist
	    VirtualParticlesSelect * pl = new ParticlesList();	    
	    try {  // try to load data
		pl->loadFile(filename,my_psv);
		bool ok;
		int step = table_list->text( i, 2 ).toInt( &ok);
		if (!ok) step=1;
		pl->step_part  = step;
                pl->defaultIndexTab();
		pl->col        = QColor( table_list->text( i, 3 ));
		pl->is_visible = ((QCheckBox*)table_list->cellWidget( i,0))->isChecked();
		ParticlesSelect * ps = new ParticlesSelect();
		ps->vps = pl;
		//pr->printRange();
		my_psv->push_back(*ps);
		PRINT_D cerr << "In SetParticlesRangeForm::apply(), my_psv->size() = "
			<<  my_psv->size() << "\n";
		delete ps;
	    } //try
	    catch (int n) {
		switch (n) {
		case -1: 
		    break;
		case -2:
		    break;	
		default:
		    assert(1);
		} //switch
		QMessageBox::information( this,QString(""),QString(pl->error_message.c_str()));
		delete pl;
	    } //catch
	}
    }
}		

// ============================================================================
// SetParticlesRangeForm::setColor()                                           
// set color at the selected cell                                              
void ParticlesSelectForm::setColor()
{
  setColor( table_range->currentRow(), table_range->currentColumn() );
  table_range->setFocus();
}

// ============================================================================
// SetParticlesRangeForm::setColor()                                           
// print color at the selected cell                                            
void ParticlesSelectForm::setColor( int row, int col )
{
  PRINT_D cerr << "Row = " << row << " col = " << col << "\n";
}
// ============================================================================
// SetParticlesRangeForm::setColor()                                           
// print color at the selected cell                                            
void ParticlesSelectForm::changeListFile( int row, int col )
{
    if (col == 1) {
	QFileInfo fileinfo(table_list->text(row,col));
	QString fn = QFileDialog::getOpenFileName( fileinfo.dirPath(TRUE), QString::null,
						   this);
	if ( !fn.isEmpty() ) {
	    table_list->setText(row,1,fn);
	}
    }
}
// ============================================================================
// SetParticlesRangeForm::currentChanged()                                     
// After a click in the 4th column (color field), color dialog box is launched 
// to select a new color                                                       
void ParticlesSelectForm::currentChanged( int row, int col, QTable * table, int col_index )
{
  colorPushButton->setEnabled( col == col_index );
  if (col == col_index ) {
    QColor color = QColorDialog::getColor(
			QColor( table->text( row, col ) ),
			this, "color dialog" );
    if ( color.isValid() ) {
	QPixmap pix = table->pixmap( row, col );
	pix.fill( color );
	table->setPixmap( row, col, pix );
	table->setText( row, col, color.name() );
	table->setCurrentCell(row,col-1); // move the focus to the previous column
    }    
  }
}
// ============================================================================
// SetParticlesRangeForm::currentChangedRange()                                     
// After a click in the 4th column (color field), color dialog box is launched 
// to select a new color                                                       
void ParticlesSelectForm::currentChangedRange( int row, int col )
{
    currentChanged( row,col,table_range,4);
}
// ============================================================================
// SetParticlesRangeForm::currentChangedRange()                                     
// After a click in the 3h column (color field), color dialog box is launched 
// to select a new color                                                       
void ParticlesSelectForm::currentChangedList( int row, int col )
{
    currentChanged( row,col,table_list,3);
}
// ============================================================================
// SetParticlesRangeForm::valueChanged()                                       
// proceed according to the changed value                                      
void ParticlesSelectForm::valueChanged( int row, int col )
{
  bool special_first_line=TRUE;
  if (col >0  && col <4  && special_first_line) { 
    bool ok;
    int value=table_range->text( row, col ).toInt( &ok );
    if (ok) {
      if ( col == 3 ) {  // step part
	table_range->setText(row, col, QString( "%1" ).arg(value,0)); // keep the value
      }
      else {
	if ( value > nbody) {
	  table_range->setText(row, col, QString( "%1" ).arg(value,0) + ">nbody!");
	}
	else {
	  if (col == 1) { // first_part
	    int other=table_range->text( row, col+1 ).toInt( &ok );
	    bool o_empty = table_range->text( row, col+1 ).isEmpty(); // is other empty ?

	    if ( ! o_empty && value > other)  { // wrong first > last
	      table_range->setText(row, col, 
                             QString( "%1" ).arg(value,0) + ">last_part!");
	    } else {
	      table_range->setText(row, col, 
                             QString( "%1" ).arg(value,0)); // keep the value
	    }
	
	  } else {        // last_part [col=2]
	    PRINT_D cerr << "current col =====> " << col << "\n";
	    int other=table_range->text( row, col-1 ).toInt( &ok );
	    bool o_empty = table_range->text( row, col-1 ).isEmpty(); // is other empty ?

	    if (! o_empty && value < other)  { // wrong last > other
	      table_range->setText(row, col, QString( "%1" ).arg(value,0) +
              "<first_part!");
	    } else {
	      table_range->setText(row, col, 
                             QString( "%1" ).arg(value,0)); // keep the value
	    }
	  }
	}
      }
    } 
    else {
      if ( !table_range->text( row, col ).isEmpty() ) {
	table_range->setText( row, col, table_range->text( row, col ) + "?" );
      }
    }
  }
}

// ============================================================================
// SetParticlesRangeForm::fillForm()                                           
// fill the form according to the data                                         
void ParticlesSelectForm::fillForm()
{
    // proceed table range
    fillFormRange();
    // proceed table list
    fillFormList();
}

// ============================================================================
// SetParticlesRangeForm::fillForm()                                           
// fill the form according to the data                                         
void ParticlesSelectForm::fillFormList()
{
  // create a pixmap to store the color
  QRect rect = table_list->cellRect( 0, 3 );
  QPixmap pix( rect.width(), rect.height() );

  // find first object vector index of type LIST
  unsigned int first_index=0;
  bool ok=false;
  for (unsigned int i=0; i<my_psv->size();i++) {
      VirtualParticlesSelect * vps = (*my_psv)[i].vps;
      first_index++;
      PRINT_D std::cerr << "in fillList, first_index =  " << first_index << " type ="<< vps->v_type <<"\n";
      if (vps->v_type==2) { // find the first index ot type LIST
	  ok=true;
	  break;
      }
  }  
  if (ok)  first_index--;   
  // fill the table_list accordung to [prv] variable
  for (unsigned int i=first_index, ii=0; i<MAX_OBJECTS+first_index; i++,ii++) {
      //std::cerr << "in fillList, first_index =  " << first_index << "\n";
      if (i < my_psv->size() && ok ) {
	  VirtualParticlesSelect * vps = (*my_psv)[i].vps;	  
	  // draw check box
	  QCheckBox * checkbox = new QCheckBox(this);
	  checkbox->setChecked(vps->is_visible);
	  table_list->setCellWidget(ii,0, checkbox);
	  // draw color and set the colorname
	  QColor color = vps->col;
	  pix.fill( color );
	  table_list->setPixmap( ii, 3, pix );   // colorize
	  table_list->setText(ii,3,vps->col.name());  // set the color name
	  // fill list file
	  table_list->setText(ii,1,QString(vps->list_file.c_str()));
	  // step part
	  table_list->setText(ii,2,QString( "%1" ).arg( vps->step_part,0));
      }
      else  {    // Fill up the remaining cells
	  QCheckBox * checkbox = new QCheckBox(this);
	  checkbox->setChecked(FALSE);
	  table_list->setCellWidget(ii,0, checkbox);
	  
	  // print range
	  table_list->setText(ii,1,"");
	  table_list->setText(ii,2,"");
	  
	  // draw black color
	  QColor color = white;
	  pix.fill( color );
	  table_list->setPixmap( ii, 3, pix );   // colorize
	  table_list->setText(ii,3,color.name());  // set the color name	  
      }	
  }
  
}
// ============================================================================
// SetParticlesRangeForm::fillForm()                                           
// fill the form according to the data                                         
void ParticlesSelectForm::fillFormRange()
{
  // create a pixmap to store the color
  QRect rect = table_range->cellRect( 0, 4 );
  QPixmap pix( rect.width(), rect.height() );

  // find first object vector index of type LIST
  unsigned int first_index=0;
  for (unsigned int i=0; i<my_psv->size();i++) {
      VirtualParticlesSelect * vps = (*my_psv)[i].vps;
      first_index++;
      if (vps->v_type==2) { // find the first index ot type LIST
	  break;
      }
  }
  
  // fill the table_range accordung to [prv] variable
  for (unsigned int i=0; i< MAX_OBJECTS; i++) {
      bool sup=true;
      bool bad=true;
      if (i < my_psv->size()) {
	  sup = false;
	  VirtualParticlesSelect * vps = (*my_psv)[i].vps;
	  if (vps->v_type==1) { // must be RANGE type
	      bad=false;
	      // draw check box
	      QCheckBox * checkbox = new QCheckBox(this);
	      checkbox->setChecked(vps->is_visible);
	      table_range->setCellWidget(i,0, checkbox);
	      
	      if (i == 0 ) { // the first row must not be editable_range to preserve 
		  // maximum particles range
		  table_range->setItem(i,1,new QTableItem( table_range, QTableItem::Never,
							   QString( "%1" ).arg( vps->first_part, 0)  ) );
		  table_range->setItem(i,2,new QTableItem( table_range, QTableItem::Never,
							   QString( "%1" ).arg( vps->last_part, 0)  ) );
		  table_range->setItem(i,3,new QTableItem( table_range, QTableItem::WhenCurrent,
							   QString( "%1" ).arg( vps->step_part, 0)  ) );
	      }	
	      else {
		  // print range
		  table_range->setText(i,1,QString( "%1" ).arg( vps->first_part, 0));
		  table_range->setText(i,2,QString( "%1" ).arg( vps->last_part, 0));
		  table_range->setText(i,3,QString( "%1" ).arg( vps->step_part, 0));
	      }		
	      // draw color and set the colorname
	      QColor color = vps->col;
	      pix.fill( color );
	      table_range->setPixmap( i, 4, pix );   // colorize
	      table_range->setText(i,4,vps->col.name());  // set the color name
	  } else {
	      bad=true;
	  }
	  
      }
      if (sup || bad)  {    // Fill up the remaining cells
	  QCheckBox * checkbox = new QCheckBox(this);
	  checkbox->setChecked(FALSE);
	  table_range->setCellWidget(i,0, checkbox);
	  
	  // print range
	  table_range->setText(i,1,"");
	  table_range->setText(i,2,"");
	  table_range->setText(i,3,"");
	  
	  // draw black color
	  QColor color = white;
	  pix.fill( color );
	  table_range->setPixmap( i, 4, pix );   // colorize
	  table_range->setText(i,4,color.name());  // set the color name
      }
  }    
}
// ============================================================================
// SetParticlesRangeForm::updateData()                                         
// fill the form according to the NEW data                                     
void ParticlesSelectForm::updateData( ParticlesSelectVector * psv, const int n_body, const float *p_pos )
{
  // save nbody for later processing
  nbody = n_body;
  // save pos for later processing
  pos   = p_pos;
  // save prv for later processing
  my_psv = psv; 
  
  fillForm();
}
// ============================================================================


