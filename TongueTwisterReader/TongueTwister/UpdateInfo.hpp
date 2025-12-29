#ifndef UpdateInfo_hpp
#define UpdateInfo_hpp

#include <map>
#include <string>
class Parameter;



struct AcceptTries {

            AcceptTries(void) { numTries = 0; numAccepts = 0; }
    int     numTries;
    int     numAccepts;
};

class UpdateInfo {

    public:
        static UpdateInfo&                                  updateInfo(void) {
                                                                static UpdateInfo uiPool;
                                                                return uiPool;
                                                                }
        void                                                accept(std::pair<Parameter*,int>& key);
        void                                                print(void);
        void                                                reject(std::pair<Parameter*,int>& key);
        AcceptTries*                                        getUpdateInfo(std::pair<Parameter*,int>& key);
        
    private:
                                                            UpdateInfo(void);
                                                           ~UpdateInfo(void);
        std::map<std::pair<Parameter*,int>,AcceptTries*>    info;
};

#endif
