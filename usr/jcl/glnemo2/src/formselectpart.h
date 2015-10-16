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
/**
	@author Jean-Charles Lambert <jean-charles.lambert@lam.fr>
 */
#ifndef GLNEMOFORMSELECTPART_H
#define GLNEMOFORMSELECTPART_H
#include "ui_formselectpart.h"
#include "snapshotinterface.h"
namespace glnemo {


class FormSelectPart: public QDialog {
  Q_OBJECT
  public:
    FormSelectPart(QWidget *parent = 0);

    ~FormSelectPart();
  void update(SnapshotInterface *,ComponentRangeVector *, const std::string, const bool first_snapshot=false);
  signals:
      void selectPart(const std::string, const bool, const bool);
  private slots:
    void on_all_check_clicked()   { updateSelect(); }
    void on_disk_check_clicked()  { updateSelect(); }
    void on_gas_check_clicked()   { updateSelect(); }
    void on_halo_check_clicked()  { updateSelect(); }
    void on_bulge_check_clicked() { updateSelect(); }
    void on_stars_check_clicked() { updateSelect(); }
    void on_bndry_check_clicked() { updateSelect(); }
    void on_load_vel_check_clicked()    { load_vel = form.load_vel_check->isChecked(); }
    void on_manual_range_textChanged(QString s) { updateSelect(); s="";}
    void on_button_box_accepted() {
      accept();
      emit selectPart(form.final_select->text().toStdString(), first_snapshot, load_vel);
    };
    void on_button_clear_clicked() { reset(false); form.manual_range->clear();}
    void on_button_clear_man_clicked() { form.manual_range->clear();}
  private:
    Ui::FormSelectPart form;
    ComponentRangeVector * crv;
    SnapshotInterface * current_data;
    bool load_vel;
    void reset(bool range=true);
    void updateSelect();
    bool first_snapshot;
};

}

#endif
