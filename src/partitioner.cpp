#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <map>
#include "cell.h"
#include "net.h"
#include "partitioner.h"
using namespace std;


void Partitioner::parseInput(fstream& inFile)
{
    string str;
    // Set balance factor
    inFile >> str;
    _bFactor = stod(str);

    // Set up whole circuit
    while (inFile >> str) {
        if (str == "NET") {
            string netName, cellName, tmpCellName = "";
            inFile >> netName;
            int netId = _netNum;
            _netArray.push_back(new Net(netName));
            _netName2Id[netName] = netId;
            while (inFile >> cellName) {
                if (cellName == ";") {
                    tmpCellName = "";
                    break;
                }
                else {
                    // a newly seen cell
                    if (_cellName2Id.count(cellName) == 0) {
                        int cellId = _cellNum;
                        _cellArray.push_back(new Cell(cellName, 0, cellId));
                        _cellName2Id[cellName] = cellId;
                        _cellArray[cellId]->addNet(netId);
                        _cellArray[cellId]->incPinNum();
                        _netArray[netId]->addCell(cellId);
                        ++_cellNum;
                        tmpCellName = cellName;
                    }
                    // an existed cell
                    else {
                        if (cellName != tmpCellName) {
                            assert(_cellName2Id.count(cellName) == 1);
                            int cellId = _cellName2Id[cellName];
                            _cellArray[cellId]->addNet(netId);
                            _cellArray[cellId]->incPinNum();
                            _netArray[netId]->addCell(cellId);
                            tmpCellName = cellName;
                        }
                    }
                }
            }
            ++_netNum;
        }
    }
    return;
}
void Partitioner::printblist(){
for (int i=-_cellNum;i<=_cellNum;i++){
            
            Node *aa=_bList[0][i],*bb=_bList[1][i];
            if (aa!=NULL)
                cout << "\ngain : "<<i<<" in A\n";
            while(aa!=NULL){
                
                if (aa->getPrev()!=NULL)
                    cout <<"|"<<_cellArray[aa->getPrev()->getId()]->getName();
                else 
                    cout <<"|NULL";
                cout<<"<-"; 
                cout <<_cellArray[aa->getId()]->getName();
                cout<<"->"; 
                if (aa->getNext()!=NULL)
                    cout <<_cellArray[aa->getNext()->getId()]->getName()<<"|\n";
                else 
                    cout <<"NULL|\n";
                aa=aa->getNext();
            }
            if(bb!=NULL)
            cout << "\ngain : "<<i<<" in B\n";
            while(bb!=NULL){
                if (bb->getPrev()!=NULL)
                    cout  <<"|"<<_cellArray[bb->getPrev()->getId()]->getName();
                else 
                    cout <<"|NULL";
                cout<<"<-"; 
                cout <<_cellArray[bb->getId()]->getName();
                cout<<"->"; 
                if (bb->getNext()!=NULL)
                    cout <<_cellArray[bb->getNext()->getId()]->getName()<<"|\n";
                else 
                    cout <<"NULL|\n";
                
                bb=bb->getNext();
            }
    }
}
void Partitioner::partition()
{
    //_cutSize
    //_partSize[0] #of cell in A
    //_partSize[1] #of cell in B
    //_maxGainCell pointer to max gain cell
    //_bList[0] A's bucket list
    //_bList[1] B's bucket list

    //_cellArray[i] 
    //             _part 0->A 1->B
    //             _gain
    //             _lock


    double L_bound,U_bound;
    L_bound=(1-_bFactor)/2*_cellNum;
    U_bound=(1+_bFactor)/2*_cellNum;
    cout <<"L_bound : "<<L_bound<<endl<<"U_bound : "<<U_bound<<endl; 
    /*****************init state*****************/
    for (int i=0;i<_cellNum;i++){
        if(i<_cellNum/2){
            _cellArray[i]->setPart(0);
            _partSize[0]++;
            for (int j=0;j< _cellArray[i]->getNetList().size();j++)
                _netArray[_cellArray[i]->getNetList()[j]]->incPartCount(0);
        }
        else{
            _cellArray[i]->setPart(1);
            _partSize[1]++;
            for (int j=0;j< _cellArray[i]->getNetList().size();j++)
                _netArray[_cellArray[i]->getNetList()[j]]->incPartCount(1);

        }
    }
    /*****************init state*****************/
    /*****************init gain*****************/
    int _A=0,_B=0;
    Node *tail[2];
    for (int i=0;i<_cellNum;i++)
        if(_cellArray[i]->getPinNum()>_maxPinNum)
            _maxPinNum=_cellArray[i]->getPinNum();
    int max_gain_value=_maxPinNum;
    for (int i=-_maxPinNum;i<=_maxPinNum;i++){
        _bList[0].insert(pair<int, Node*>(i,NULL));
        _bList[1].insert(pair<int, Node*>(i,NULL));
    }
    for (int i=0;i<_cellNum;i++){
        if(_cellArray[i]->getPart()==1)
            for (int j=0;j<_cellArray[i]->getNetList().size();j++){
                if(_netArray[_cellArray[i]->getNetList()[j]]->getPartCount(1)==1)//from size ==1  -> gain+1
                    _cellArray[i]->incGain();
                else if(_netArray[_cellArray[i]->getNetList()[j]]->getPartCount(0)==0)//to size == 0 -> gain-1
                    _cellArray[i]->decGain();
            }
        else
            for (int j=0;j<_cellArray[i]->getNetList().size();j++){
                if(_netArray[_cellArray[i]->getNetList()[j]]->getPartCount(0)==1)//from size ==1 -> gain+1
                    _cellArray[i]->incGain();
                else if(_netArray[_cellArray[i]->getNetList()[j]]->getPartCount(1)==0)//to size == 0 -> gain-1
                    _cellArray[i]->decGain();
            }
        Node *n=new Node(i);
        _cellArray[i]->setNode(n);
        if(max_gain_value<_cellArray[i]->getGain())
            max_gain_value=_cellArray[i]->getGain();
        int p=_cellArray[i]->getPart();
        if(_bList[p][_cellArray[i]->getGain()]==NULL){
            _bList[p][_cellArray[i]->getGain()]=n;
        }
        else{
            _bList[p][_cellArray[i]->getGain()]->setPrev(n);
            n->setNext(_bList[p][_cellArray[i]->getGain()]);
            _bList[p][_cellArray[i]->getGain()]=n;
        }
    }
    /*output gain
    for (int i=-_cellNum;i<=_cellNum;i++){
        cout << "\ngain : "<<i<<" in A";
        Node *aa=_bList[0][i],*bb=_bList[1][i];
        while(aa!=NULL){
            cout <<aa->getId()<<" ";
            aa=aa->getNext();
        }
        cout << "\ngain : "<<i<<" in B";
        while(bb!=NULL){
            cout <<bb->getId()<<" ";
            bb=bb->getNext();
        }
    }
    */
    /*****************init gain*****************/
    /*****************init cutsize*****************/
    for(int i=0;i<_netArray.size();i++)
        if(_netArray[i]->getPartCount(0)>0&&_netArray[i]->getPartCount(1)>0)
            _cutSize++;
    cout <<"\ncutsize : "<<_cutSize<<endl;
    /*****************init cutsize*****************/
    /*****************process start*****************/

    //select a node with max gain and more balance
    
    bool flag=true;
    cout <<"======step start======"<<endl;
    while(flag){
        //printblist(); 
        _maxGainCell=NULL;

        //get max gainpointer
        while(true){
            if(_bList[0][max_gain_value]!=NULL)
                _maxGainCell=_bList[0][max_gain_value];
            while(_maxGainCell!=NULL){ 
                if (L_bound < _partSize[_cellArray[_maxGainCell->getId()]->getPart()] && _partSize[_cellArray[_maxGainCell->getId()]->getPart()] < U_bound)
                    break;
                
                _maxGainCell=_maxGainCell->getNext();
            }

            if(_maxGainCell==NULL)
                _maxGainCell=_bList[1][max_gain_value];
            while(_maxGainCell!=NULL){
                if (L_bound < _partSize[_cellArray[_maxGainCell->getId()]->getPart()] && _partSize[_cellArray[_maxGainCell->getId()]->getPart()] < U_bound)
                    break;
                
                _maxGainCell=_maxGainCell->getNext();
            }
            if(_maxGainCell!=NULL || max_gain_value< -_maxPinNum)
                break;
            max_gain_value--;
        }

        //if gainpointer not fount(means that is the end of the iteration)
        if(_maxGainCell==NULL){
             cout <<"======iteration end======\n\n";
             cout <<"move his : ";
             for (int i=0;i<_moveStack.size();i++)
                cout << _cellArray[_moveStack[i]]->getName()<<" ";
            cout <<endl<<"_moveNum : "<<_moveNum<<endl;
            cout <<"_bestMoveNum : "<<_bestMoveNum<<endl;
            cout <<"bestMove : ";
            for (int i=0;i<=_bestMoveNum;i++)
                cout << _cellArray[_moveStack[i]]->getName()<<" ";
            break;
        }

        //remove node from bucket
        if (_maxGainCell->getPrev()!=NULL)
            _maxGainCell->getPrev()->setNext(_maxGainCell->getNext());
        else
            _bList[_cellArray[_maxGainCell->getId()]->getPart()][_cellArray[_maxGainCell->getId()]->getGain()]=_maxGainCell->getNext();
            
        if (_maxGainCell->getNext()!=NULL)
            _maxGainCell->getNext()->setPrev(_maxGainCell->getPrev());
        
        //lock the switched node
        _cellArray[_maxGainCell->getId()]->lock();

        cout <<"\nselect node : "<<_cellArray[_maxGainCell->getId()]->getName()<<endl
             <<"gain : "<<_cellArray[_maxGainCell->getId()]->getGain()<<endl;

        update_gain(_cellArray[_maxGainCell->getId()],max_gain_value);

        
        
        _accGain+=_cellArray[_maxGainCell->getId()]->getGain();
        cout <<"_accGain : "<<_accGain<<endl;
        if (_accGain>_maxAccGain){
            _maxAccGain=_accGain;
            _bestMoveNum=_moveNum;
        }

        //cout << _maxGainCell->getId()<<endl;
        _moveNum++;
        _moveStack.push_back(_maxGainCell->getId());
    }
    /*****************process end*****************/
    //reportNet() ;
    //reportCell(); 
}
void Partitioner::update_gain(Cell* target,int &max_gain_value){
    bool part=target->getPart();
    vector <int> ori_gain;
    max_gain_value=_maxPinNum;
    for (int i=0;i<_cellNum;i++)
        ori_gain.push_back(_cellArray[i]->getGain());
    for (int i=0;i<target->getNetList().size();i++){
        if(_netArray[target->getNetList()[i]]->getPartCount(part)==1){
            for(int j=0;j<_netArray[target->getNetList()[i]]->getCellList().size();j++){
                if(_cellArray[_netArray[target->getNetList()[i]]->getCellList()[j]]->getPart()==!part && _cellArray[_netArray[target->getNetList()[i]]->getCellList()[j]]->getName()!=target->getName()){
                    _cellArray[_netArray[target->getNetList()[i]]->getCellList()[j]]->decGain();
                }
            }
        }
        if(_netArray[target->getNetList()[i]]->getPartCount(part)==2){
            for(int j=0;j<_netArray[target->getNetList()[i]]->getCellList().size();j++){
                if(_cellArray[_netArray[target->getNetList()[i]]->getCellList()[j]]->getPart()==part && _cellArray[_netArray[target->getNetList()[i]]->getCellList()[j]]->getName()!=target->getName())
                    _cellArray[_netArray[target->getNetList()[i]]->getCellList()[j]]->incGain();
            }
        }
        if(_netArray[target->getNetList()[i]]->getPartCount(!part)==0){
            for(int j=0;j<_netArray[target->getNetList()[i]]->getCellList().size();j++){
                if(_cellArray[_netArray[target->getNetList()[i]]->getCellList()[j]]->getPart()==part && _cellArray[_netArray[target->getNetList()[i]]->getCellList()[j]]->getName()!=target->getName())
                    _cellArray[_netArray[target->getNetList()[i]]->getCellList()[j]]->incGain();
            }
        }
        if(_netArray[target->getNetList()[i]]->getPartCount(!part)==1){
            for(int j=0;j<_netArray[target->getNetList()[i]]->getCellList().size();j++){
                if(_cellArray[_netArray[target->getNetList()[i]]->getCellList()[j]]->getPart()==!part && _cellArray[_netArray[target->getNetList()[i]]->getCellList()[j]]->getName()!=target->getName())
                    _cellArray[_netArray[target->getNetList()[i]]->getCellList()[j]]->decGain();
            }
        }
    }
    for (int i=0;i<_cellNum; i++){
        if(_cellArray[i]->getGain()>max_gain_value)
            max_gain_value=_cellArray[i]->getGain();
        if(ori_gain[i]!=_cellArray[i]->getGain() && !_cellArray[i]->getLock() ){
            
            if(_cellArray[i]->getNode()->getNext()!=NULL)
                _cellArray[i]->getNode()->getNext()->setPrev(_cellArray[i]->getNode()->getPrev());
            if(_cellArray[i]->getNode()->getPrev()!=NULL)
                _cellArray[i]->getNode()->getPrev()->setNext(_cellArray[i]->getNode()->getNext());
            if(_cellArray[i]->getNode()->getPrev()==NULL){
                _bList[_cellArray[i]->getPart()][ori_gain[i]]=_cellArray[i]->getNode()->getNext();
            }
            if(_bList[_cellArray[i]->getPart()][_cellArray[i]->getGain()]!=NULL){
                _bList[_cellArray[i]->getPart()][_cellArray[i]->getGain()]->setPrev(_cellArray[i]->getNode());
                _cellArray[i]->getNode()->setNext(_bList[_cellArray[i]->getPart()][_cellArray[i]->getGain()]);
            }
            _bList[_cellArray[i]->getPart()][_cellArray[i]->getGain()]=_cellArray[i]->getNode();
            _bList[_cellArray[i]->getPart()][_cellArray[i]->getGain()]->setPrev(NULL);
        }


    }
   //printblist();
}
void Partitioner::printSummary() const
{
    cout << endl;
    cout << "==================== Summary ====================" << endl;
    cout << " Cutsize: " << _cutSize << endl;
    cout << " Total cell number: " << _cellNum << endl;
    cout << " Total net number:  " << _netNum << endl;
    cout << " Cell Number of partition A: " << _partSize[0] << endl;
    cout << " Cell Number of partition B: " << _partSize[1] << endl;
    cout << "=================================================" << endl;
    cout << endl;
    return;
}

void Partitioner::reportNet() const
{
    cout << "Number of nets: " << _netNum << endl;
    for (size_t i = 0, end_i = _netArray.size(); i < end_i; ++i) {
        cout << setw(8) << _netArray[i]->getName() << ": ";
        vector<int> cellList = _netArray[i]->getCellList();
        for (size_t j = 0, end_j = cellList.size(); j < end_j; ++j) {
            cout << setw(8) << _cellArray[cellList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::reportCell() const
{
    cout << "Number of cells: " << _cellNum << endl;
    for (size_t i = 0, end_i = _cellArray.size(); i < end_i; ++i) {
        cout << setw(8) << _cellArray[i]->getName() << ": ";
        vector<int> netList = _cellArray[i]->getNetList();
        for (size_t j = 0, end_j = netList.size(); j < end_j; ++j) {
            cout << setw(8) << _netArray[netList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::writeResult(fstream& outFile)
{
    stringstream buff;
    buff << _cutSize;
    outFile << "Cutsize = " << buff.str() << '\n';
    buff.str("");
    buff << _partSize[0];
    outFile << "G1 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 0) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    buff.str("");
    buff << _partSize[1];
    outFile << "G2 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 1) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    return;
}

void Partitioner::clear()
{
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        delete _cellArray[i];
    }
    for (size_t i = 0, end = _netArray.size(); i < end; ++i) {
        delete _netArray[i];
    }
    return;
}
