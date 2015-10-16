// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2015                                  
// e-mail:   Jean-Charles.Lambert@lam.fr                                      
// address:  Centre de donneeS Astrophysique de Marseille (CeSAM)              
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#include "formselectpart.h"
#include <sstream>
namespace glnemo {

// ============================================================================
// Constructor                                                                 
FormSelectPart::FormSelectPart(QWidget *parent)
{
  if (parent) {;}  // remove compiler warning
  form.setupUi(this);
  load_vel = false;
}

// ============================================================================
// Destructor                                                                  
FormSelectPart::~FormSelectPart()
{
}
// ============================================================================
// update                                                                      
void FormSelectPart::update(SnapshotInterface * _si,
                            ComponentRangeVector * _crv, 
                            const std::string previous_sel, const bool first)
{
  current_data = _si;
  crv = _crv;
  std::cerr << " FormSelectPart::update ----->  \n";
  ComponentRange::list(crv);
  first_snapshot = first;
  // update nbody,time and data type
  std::ostringstream stm1,stm2;
  stm1 << current_data->nbody_first; // nbody
  form.nbody->setText((stm1.str()).c_str());
  stm2 << current_data->time_first;  // time
  form.time->setText((stm2.str()).c_str());
  form.data_type->setText(QString((current_data->getInterfaceType()).c_str())); // data type
  reset(); // reset fields
  
  bool exist_sel=false;
  if (! previous_sel.empty()) exist_sel=true;
  form.manual_range->setText(QString(previous_sel.c_str())); // previous selected range
  // loop on all possible components
  for (std::vector<ComponentRange>::iterator icrv=crv->begin(); icrv<crv->end();icrv++) {
    if (icrv->type == "all" ) {
      if (!exist_sel && crv->size()==1) form.all_check->setChecked(true);
      form.all_range->setAlignment(Qt::AlignLeft);
      form.all_range->setText(QString((icrv->range).c_str()));
    }
    if (icrv->type == "disk" ) {
      if (!exist_sel) form.disk_check->setChecked(true);
      form.disk_range->setAlignment(Qt::AlignLeft);
      form.disk_range->setText(QString((icrv->range).c_str()));
     }
    if (icrv->type == "halo" ) {
      if (!exist_sel) form.halo_check->setChecked(true);
      form.halo_range->setAlignment(Qt::AlignLeft);
      form.halo_range->setText(QString((icrv->range).c_str()));
    }
    if (icrv->type == "gas" ) {
      if (!exist_sel) form.gas_check->setChecked(true);
      form.gas_range->setAlignment(Qt::AlignLeft);
      form.gas_range->setText(QString((icrv->range).c_str()));
    }
    if (icrv->type == "bulge" ) {
      if (!exist_sel) form.bulge_check->setChecked(true);
      form.bulge_range->setAlignment(Qt::AlignLeft);
      form.bulge_range->setText(QString((icrv->range).c_str()));
    }
    if (icrv->type == "stars" ) {
      if (!exist_sel) form.stars_check->setChecked(true);
      form.stars_range->setAlignment(Qt::AlignLeft);
      form.stars_range->setText(QString((icrv->range).c_str()));
    }
    if (icrv->type == "bndry" ) {
      if (!exist_sel) form.bndry_check->setChecked(true);
      form.bndry_range->setAlignment(Qt::AlignLeft);
      form.bndry_range->setText(QString((icrv->range).c_str()));
    }
  }
  updateSelect();
}

// ============================================================================
// reset                                                                       
// reset all the fiels of the form                                             
void FormSelectPart::reset(bool range)
{
  form.all_check->setChecked(false);
  if (range) {
    form.all_range->setAlignment(Qt::AlignRight);
    form.all_range->setText("Unknown");
  }

  form.disk_check->setChecked(false);
  if (range) {
    form.disk_range->setAlignment(Qt::AlignRight);
    form.disk_range->setText("Unknown");
  }

  form.halo_check->setChecked(false);
  if (range) {
    form.halo_range->setAlignment(Qt::AlignRight);
    form.halo_range->setText("Unknown");
  }

  form.gas_check->setChecked(false);
  if (range) {
    form.gas_range->setAlignment(Qt::AlignRight);
    form.gas_range->setText("Unknown");
  }

  form.bulge_check->setChecked(false);
  if (range) {
    form.bulge_range->setAlignment(Qt::AlignRight);
    form.bulge_range->setText("Unknown");
  }

  form.stars_check->setChecked(false);
  if (range) {
    form.stars_range->setAlignment(Qt::AlignRight);
    form.stars_range->setText("Unknown");
  }

  form.bndry_check->setChecked(false);
  if (range) {
    form.bndry_range->setAlignment(Qt::AlignRight);
    form.bndry_range->setText("Unknown");
  }
  form.final_select->clear();
}
// ============================================================================
// updateSelect                                                                
void FormSelectPart::updateSelect()
{
  bool coma=true;
  QString tcoma="";
  form.final_select->setText(form.manual_range->text());
  if (form.final_select->text().isEmpty()) coma=false;
  if (coma) tcoma=",";
  if (form.all_check->isChecked()) {
    form.final_select->insert(tcoma+QString("all"));
    coma=true;
  }
  if (coma) tcoma=",";
  if (form.halo_check->isChecked()) {
    form.final_select->insert(tcoma+QString("halo"));
    coma=true;
  }  
  if (coma) tcoma=",";
  if (form.gas_check->isChecked()) {
    form.final_select->insert(tcoma+QString("gas"));
    coma=true;
  }
  if (coma) tcoma=",";
  if (form.disk_check->isChecked()) {
    form.final_select->insert(tcoma+QString("disk"));
    coma=true;
  }
  if (coma) tcoma=",";
  if (form.stars_check->isChecked()) {
    form.final_select->insert(tcoma+QString("stars"));
    coma=true;
  }
  if (coma) tcoma=",";
  if (form.bulge_check->isChecked()) {
    form.final_select->insert(tcoma+QString("bulge"));
    coma=true;
  }
  if (coma) tcoma=",";
  if (form.bndry_check->isChecked()) {
    form.final_select->insert(tcoma+QString("bndry"));
    coma=true;
  }
}
}
