// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
//                                                                             
// Implementation of the form which set the particles range & color            
//                                                                             
//                                                                             
// ============================================================================

#include <qcolordialog.h>
#include <qcombobox.h>
#include <qlayout.h>
#include <qpixmap.h>
#include <qpushbutton.h>
#include <qtable.h>
#include <qcheckbox.h>

#include "set_particles_range_form.h"

#define LOCAL_DEBUG 0
#include "print_debug.h"

const unsigned int MAX_OBJECTS=30;
// ============================================================================
// constructor
SetParticlesRangeForm::SetParticlesRangeForm(const GLBox * glbox,
					     ParticlesRangeVector * prv, 
                                             const int n_body,
					     const float * p_pos,
					     QWidget *parent , const char *name,
					     bool modal, WFlags f ): 
                                             QDialog( parent, name, modal, f )
{
  setCaption( "Set particles range and color" );
  //resize( 400, 240 );
  resize( 370, 174 );
  tableButtonBox = new QVBoxLayout( this, 11, 6, "table button box layout" );

  // save nbody for later processing
  nbody = n_body;
  
  // save pos for later processing
  pos   = p_pos;

  // save prv for later processing
  my_prv = prv;

  // Create the QTable
  table = new QTable( this, "data table" );
  table->setNumCols( 5 );
  table->setNumRows(MAX_OBJECTS);
  //table->setColumnReadOnly( 1, TRUE );
  //table->setColumnReadOnly( 2, TRUE );
  //table->setColumnReadOnly( 4, TRUE );
  table->setColumnWidth( 0, 50 );
  table->setColumnWidth( 1, 80 );
  table->setColumnWidth( 2, 80 );
  table->setColumnWidth( 3, 40 );
  table->setColumnWidth( 4, 50 );
  QHeader *th = table->horizontalHeader();
  th->setLabel( 0, "Visible" );
  th->setLabel( 1, "First part" );
  th->setLabel( 2, "Last part" );
  th->setLabel( 3, "Step" );
  th->setLabel( 4, "Color" );
  tableButtonBox->addWidget( table );  // insert table in the layout

  // button box
  buttonBox = new QHBoxLayout( 0, 0, 6, "button box layout" );

  // color button
  colorPushButton = new QPushButton( this, "color button" );
  colorPushButton->setText( "&Color..." );
  colorPushButton->setEnabled( FALSE );
  buttonBox->addWidget( colorPushButton );


  QSpacerItem *spacer = new QSpacerItem( 0, 0, QSizePolicy::Expanding,
					 QSizePolicy::Minimum );
  buttonBox->addItem( spacer );

  okPushButton = new QPushButton( this, "ok button" );
  okPushButton->setText( "OK" );
  okPushButton->setDefault( TRUE );
  buttonBox->addWidget( okPushButton );

  applyPushButton = new QPushButton( this, "apply button" );
  applyPushButton->setText( "Apply" );
  applyPushButton->setDefault( TRUE );
  buttonBox->addWidget( applyPushButton );

  cancelPushButton = new QPushButton( this, "cancel button" );
  cancelPushButton->setText( "Cancel" );
  cancelPushButton->setAccel( Key_Escape );
  buttonBox->addWidget( cancelPushButton );

  tableButtonBox->addLayout( buttonBox ); // insert buttonBox in the layout

  // Connexion with the TABLE widget
  connect( table, SIGNAL( clicked(int,int,int,const QPoint&) ),
	   this, SLOT( setColor(int,int) ) );
  connect( table, SIGNAL( currentChanged(int,int) ),
	   this, SLOT( currentChanged(int,int) ) );
  connect( table, SIGNAL( valueChanged(int,int) ),
	   this, SLOT( valueChanged(int,int) ) );
  connect( colorPushButton, SIGNAL( clicked() ), this, SLOT( setColor() ) );
  connect( okPushButton, SIGNAL( clicked() ), this, SLOT( accept()) );
  connect( applyPushButton, SIGNAL( clicked() ), this, SLOT( apply()) );
  connect( cancelPushButton, SIGNAL( clicked() ), this, SLOT( reject() ) );

  // connection to update the data
    connect(this,SIGNAL(applyData(const int *, const float *, 
				const ParticlesRangeVector * )),
	  glbox,SLOT(getData(const int *, const float *, 
                     const ParticlesRangeVector *)));

  // fill FORM
  fillForm();
  
  clearWState( WState_Polished ); 
}
// ============================================================================
// destructor
SetParticlesRangeForm::~SetParticlesRangeForm()
{
}
// ============================================================================
// fill the form according to the NEW data
void SetParticlesRangeForm::updateData(ParticlesRangeVector*prv,
                                       const int n_body,const float *p_pos)
{
   // save nbody for later processing
  nbody = n_body;
  
  // save pos for later processing
  pos   = p_pos;

  // save prv for later processing
  my_prv = prv; 
  
  fillForm();
}                                       
// ============================================================================
// fill the form according to the data
void SetParticlesRangeForm::fillForm()
{

  // create a pixmap to store the color
  QRect rect = table->cellRect( 0, 4 );
  QPixmap pix( rect.width(), rect.height() );

  // fill the table accordung to [prv] variable
  for (unsigned int i=0; i< MAX_OBJECTS; i++) {
    if (i < my_prv->size()) {
      ParticlesRange pr = (* my_prv)[i];

      // draw check box
      QCheckBox * checkbox = new QCheckBox(this);
      checkbox->setChecked(pr.is_visible);
      table->setCellWidget(i,0, checkbox);

      if (i == 0 ) { // the first row must not be editable to preserve 
                     // maximum particles range
        table->setItem(i,1,new QTableItem( table, QTableItem::Never,
                               QString( "%1" ).arg( pr.first_part, 0)  ) );
        table->setItem(i,2,new QTableItem( table, QTableItem::Never,
                               QString( "%1" ).arg( pr.last_part, 0)  ) );
        table->setItem(i,3,new QTableItem( table, QTableItem::WhenCurrent,
                               QString( "%1" ).arg( pr.step_part, 0)  ) );
      }
      else {
        // print range
        table->setText(i,1,QString( "%1" ).arg( pr.first_part, 0));
        table->setText(i,2,QString( "%1" ).arg( pr.last_part, 0));
        table->setText(i,3,QString( "%1" ).arg( pr.step_part, 0));
      }
      // draw color and set the colorname
      QColor color = pr.col;
      pix.fill( color );
      table->setPixmap( i, 4, pix );   // colorize
      table->setText(i,4,pr.col.name());  // set the color name
    }
    else {    // Fill up the remaining cells
      QCheckBox * checkbox = new QCheckBox(this);
      checkbox->setChecked(FALSE);
      table->setCellWidget(i,0, checkbox);
       
       // print range
      table->setText(i,1,"");
      table->setText(i,2,"");
      table->setText(i,3,"");
     
      // draw black color
      QColor color = white;
      pix.fill( color );
      table->setPixmap( i, 4, pix );   // colorize
      table->setText(i,4,color.name());  // set the color name
    }
  }
}
// ============================================================================
// if button color is clicked
void SetParticlesRangeForm::setColor()
{
  setColor( table->currentRow(), table->currentColumn() );
  table->setFocus();
}
// ============================================================================
// 
void SetParticlesRangeForm::setColor(int row, int col )
{
  PRINT_D cerr << "Row = " << row << " col = " << col << "\n";
}
// ============================================================================
// After a click in the 4th column (color field), color dialog box is launched
// to select a new color
void SetParticlesRangeForm::currentChanged( int row, int col )
{
  colorPushButton->setEnabled( col == 4 );
  if (col == 4 ) {
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
// 
void SetParticlesRangeForm::valueChanged( int row, int col )
{
  bool special_first_line=TRUE;
#if 0  
  if (row ==0 && col != 0) {
    special_first_line = FALSE; // Does NOT grant first line if visible NOT selected
  }
#endif
  if (col >0  && col <4  && special_first_line) { 
    bool ok;
    int value=table->text( row, col ).toInt( &ok );
    if (ok) {
      if ( col == 3 ) {  // step part
	table->setText(row, col, QString( "%1" ).arg(value,0)); // keep the value
      }
      else {
	if ( value > nbody) {
	  table->setText(row, col, QString( "%1" ).arg(value,0) + ">nbody!");
	}
	else {
	  if (col == 1) { // first_part
	    int other=table->text( row, col+1 ).toInt( &ok );
	    bool o_empty = table->text( row, col+1 ).isEmpty(); // is other empty ?

	    if ( ! o_empty && value > other)  { // wrong first > last
	      table->setText(row, col, 
                             QString( "%1" ).arg(value,0) + ">last_part!");
	    } else {
	      table->setText(row, col, 
                             QString( "%1" ).arg(value,0)); // keep the value
	    }
	
	  } else {        // last_part [col=2]
	    PRINT_D cerr << "current col =====> " << col << "\n";
	    int other=table->text( row, col-1 ).toInt( &ok );
	    bool o_empty = table->text( row, col-1 ).isEmpty(); // is other empty ?

	    if (! o_empty && value < other)  { // wrong last > other
	      table->setText(row, col, QString( "%1" ).arg(value,0) +
              "<first_part!");
	    } else {
	      table->setText(row, col, 
                             QString( "%1" ).arg(value,0)); // keep the value
	    }
	  }
	}
      }
    } 
    else {
      if ( !table->text( row, col ).isEmpty() ) {
	table->setText( row, col, table->text( row, col ) + "?" );
      }
    }
  }
}
// ============================================================================
// 
void SetParticlesRangeForm::apply()
{
  my_prv->clear(); // Clear the vector
  ParticlesRange::nb_select=0;

  for (unsigned int i=0; i < MAX_OBJECTS; i++) {
    bool ok1,ok2,commit=FALSE;
    int first = table->text( i, 1 ).toInt( &ok1 );
    int last  = table->text( i, 2 ).toInt( &ok2 );

    int step;
    if ( ok1 && ok2 ) {
      commit = TRUE;
      step= table->text( i, 3 ).toInt( &ok2 );
      if (!ok2) step=1;    // force step to 1
    }
    if (commit) {
      ParticlesRange * pr = new ParticlesRange();
      pr->first_part = first;
      pr->last_part  = last;
      pr->step_part  = step;
      pr->col        = QColor( table->text( i, 4 ));
      pr->is_visible = ((QCheckBox*)table->cellWidget( i,0))->isChecked();
      pr->printRange();
      my_prv->push_back(*pr);
      PRINT_D cerr << "In globwin, my_prv->size() = " << my_prv->size() << "\n";
      delete pr;
    }
  }
  emit applyData(&nbody,pos,my_prv);  // send signal to GLBox to update the data
}
// ============================================================================
// 
void SetParticlesRangeForm::accept()
{
  apply();
  QDialog::accept();
}
// ============================================================================
// 
