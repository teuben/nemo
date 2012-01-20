// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2012                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
/**
	@author Jean-Charles Lambert <Jean-Charles.Lambert@oamp.fr>
 */

#ifndef GLNEMOSNAPSHOTINTERFACE_H
#define GLNEMOSNAPSHOTINTERFACE_H
#include <iostream>
#include <assert.h>
#include <cmath>
#include <QFile>
#include "particlesobject.h"
#include "particlesdata.h"
#include "componentrange.h"
#include "globaloptions.h"
#include <QObject>

namespace glnemo {
//class ParticlesData;
  class CSelectTime;
  typedef std::vector<CSelectTime> CSelectTimeVector;

  //                             
  // Class to store selected time
  //                             
  class CSelectTime {
  public:
    CSelectTime(const float _i, const float _s, 
		const float _o, const float _l):inf(_i),sup(_s),offset(_o),lastt(_l) {
    };
    const CSelectTime & operator=(const CSelectTime& m) {
      inf=m.inf;
      sup=m.sup;
      lastt=m.lastt;
      offset=m.offset;
      return *this;
    }
    float inf,sup,offset,lastt;
  };
class SnapshotInterface: public QObject
{
  Q_OBJECT
  public:
    SnapshotInterface()
    {
      obj=NULL; pos=NULL; vel=NULL;
      part_data=NULL;
      end_of_data=false;
      first=true;
      crv.clear();
      stv.clear();
      parseSelectTime();
    }
    virtual ~SnapshotInterface() {};
    // ---------------------------------------------------
    // Pure Virtual functions, *** MUST be implemented ***
    // ---------------------------------------------------
    virtual SnapshotInterface * newObject(const std::string _filename, const int x=0) = 0 ;
    virtual bool isValidData() = 0;
    virtual ComponentRangeVector * getSnapshotRange() = 0;
    virtual int initLoading(GlobalOptions * so) = 0;
    virtual int nextFrame(const int * index_tab, const int nsel)= 0;
           // index_tab = array of selected indexes (size max=part_data->nbody)
           // nsel      = #particles selected in index_tab                     
           // particles not selected must have the value '-1'                  
    virtual int close() = 0;
    virtual QString endOfDataMessage() = 0;
    // simple Virtual functions
	virtual bool isEndOfData() const { return end_of_data;}
	virtual void setPort(const int x=4000) { port=x;}
	virtual void setSelectPart(const std::string _sel) { select_part = _sel;}
	virtual std::string getSelectPart() { return select_part; }
    // normal functions        
	void setFileName(std::string _f) { filename = _f;}
	std::string getFileName() const { return filename;}
	std::string getInterfaceType() { return interface_type;}
	bool isFileExist() { return QFile::exists(QString(filename.c_str()));}

    ParticlesData * part_data;
    int nbody_first;
    float time_first;
    ComponentRangeVector crv_first;
signals:
    void stringStatus(const QString);
    void intStatus(const int);
    
public slots:
    void slotStringStatus(const QString status) {
      emit stringStatus(status);
    }
    void slotIntStatus(const int status) {
      emit intStatus(status);
    }
protected:
    SnapshotInterface * obj;
    std::string filename;
    std::string interface_type;
    mutable bool end_of_data;
    std::string select_part, select_time;
    int port;
    bool keep_all, load_vel;
    ComponentRangeVector crv;
    float * pos, *vel;
    bool first;
    CSelectTimeVector stv;
    void parseSelectTime();
    std::string parseString(std::string & next_string);
    inline bool diffTime(const float t,const float fuzz=0.000001) {
      return ((fabs(t)<fuzz)?true:false);
    }
    void getRangeTime(std::string);
    bool checkRangeTime(const float);

};

}
Q_DECLARE_INTERFACE(glnemo::SnapshotInterface,
                    "Glnemo.SnapshotInterface/1.0");
#endif
