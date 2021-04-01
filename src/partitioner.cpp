#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
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

void Partitioner::printblist() {

	cout << "===============\n";
	for (int i = -_cellNum; i <= _cellNum; i++) {
		if (i == 2) {
			Node* aa = _bList[0][i], * bb = _bList[1][i];
			if (aa != NULL)
				cout << "\ngain : " << i << " in A\n";
			while (aa != NULL) {
				if (aa->getPrev() != NULL)
					cout << "|" << _cellArray[aa->getPrev()->getId()]->getName();
				else
					cout << "|NULL";
				cout << "<-";
				cout << _cellArray[aa->getId()]->getName();
				cout << "->";
				if (aa->getNext() != NULL)
					cout << _cellArray[aa->getNext()->getId()]->getName() << "|\n";
				else
					cout << "NULL|\n";
				aa = aa->getNext();
			}
			if (bb != NULL)
				cout << "\ngain : " << i << " in B\n";
			while (bb != NULL) {
				if (bb->getPrev() != NULL)
					cout << "|" << _cellArray[bb->getPrev()->getId()]->getName();
				else
					cout << "|NULL";
				cout << "<-";
				cout << _cellArray[bb->getId()]->getName();
				cout << "->";
				if (bb->getNext() != NULL)
					cout << _cellArray[bb->getNext()->getId()]->getName() << "|\n";
				else
					cout << "NULL|\n";

				bb = bb->getNext();
			}
		}
	}
}

void Partitioner::init_gain(int &max_gain_value) {

	for (int i = 0; i < _cellNum; i++) {
		_cellArray[i]->setGain(0);
		if (_cellArray[i]->getPart() == 1)
			for (int j = 0; j < _cellArray[i]->getNetList().size(); j++) {
				if (_netArray[_cellArray[i]->getNetList()[j]]->getPartCount(1) == 1)//from size ==1  -> gain+1
					_cellArray[i]->incGain();
				 if (_netArray[_cellArray[i]->getNetList()[j]]->getPartCount(0) == 0)//to size == 0 -> gain-1
					_cellArray[i]->decGain();
			}
		else
			for (int j = 0; j < _cellArray[i]->getNetList().size(); j++) {
				if (_netArray[_cellArray[i]->getNetList()[j]]->getPartCount(0) == 1)//from size ==1 -> gain+1
					_cellArray[i]->incGain();
				if (_netArray[_cellArray[i]->getNetList()[j]]->getPartCount(1) == 0)//to size == 0 -> gain-1
					_cellArray[i]->decGain();
			}
		Node* n = new Node(i);
		_cellArray[i]->setNode(n);
		if (_cellArray[i]->getGain() > max_gain_value)
			max_gain_value = _cellArray[i]->getGain();
		int p = _cellArray[i]->getPart();
		if (_bList[p][_cellArray[i]->getGain()] == NULL) {
			_bList[p][_cellArray[i]->getGain()] = n;
		}
		else {
			_bList[p][_cellArray[i]->getGain()]->setPrev(n);
			n->setNext(_bList[p][_cellArray[i]->getGain()]);
			_bList[p][_cellArray[i]->getGain()] = n;
		}
	}
}

void Partitioner::partition()
{
	double L_bound, U_bound;
	L_bound = (1 - _bFactor) / 2 * _cellNum;
	U_bound = (1 + _bFactor) / 2 * _cellNum;
	cout << "L_bound : " << L_bound << endl << "U_bound : " << U_bound << endl;

	for (int i = 0; i < _cellNum; i++) {
		if (_cellArray[i]->getPinNum() > _maxPinNum)
			_maxPinNum = _cellArray[i]->getPinNum();
	}
	for (int i = -_maxPinNum; i <= _maxPinNum; i++) {
		_bList[0].insert(pair<int, Node*>(i, NULL));
		_bList[1].insert(pair<int, Node*>(i, NULL));
	}
	/*****************init state start*****************/
	for (int i = 0; i < _cellNum; i++) {
		if (i < _cellNum / 2) {
			_cellArray[i]->setPart(0);
			_partSize[0]++;
			for (int j = 0; j < _cellArray[i]->getNetList().size(); j++)
				_netArray[_cellArray[i]->getNetList()[j]]->incPartCount(0);
		}
		else {
			_cellArray[i]->setPart(1);
			_partSize[1]++;
			for (int j = 0; j < _cellArray[i]->getNetList().size(); j++)
				_netArray[_cellArray[i]->getNetList()[j]]->incPartCount(1);

		}
	}
	/*****************init state end*****************/

	/*****************init cutsize start*****************/
	for (int i = 0; i < _netArray.size(); i++)
		if (_netArray[i]->getPartCount(0) > 0 && _netArray[i]->getPartCount(1) > 0)
			_cutSize++;
	/*****************init cutsize end*****************/

	/*****************process start*****************/
	//select a node with max gain and more balance
	while (true) {
		//cout << _iterNum << endl;
		_iterNum++;
		int max_gain_value = -_maxPinNum;
		_maxAccGain = 0;
		_accGain = 0;
		init_gain(max_gain_value);
		while (true) {
			_maxGainCell = NULL;
			while (true) {
				_maxGainCell = _bList[0][max_gain_value];
				if (L_bound <= _partSize[0] - 1 && _partSize[0] - 1 <= U_bound && _maxGainCell!=NULL)
					break;
				_maxGainCell = _bList[1][max_gain_value];
				if (L_bound <= _partSize[1] - 1 && _partSize[1] - 1 <= U_bound && _maxGainCell != NULL )
					break;
				if (max_gain_value < -_maxPinNum)
					break;
				max_gain_value--;
			}
			//if gainpointer not fount(means that is the end of the iteration)
			if (_maxGainCell == NULL) {
				if (_maxAccGain <= 0) {
					swap(_partSize[0], _partSize[1]);
					return;
				}
				if (_maxAccGain > 0) {
					for (int i = _bestMoveNum+1; i <_moveStack.size() ; i++) {
						for (int j = 0; j < _cellArray[_moveStack[i]]->getNetList().size(); j++) {
							_netArray[_cellArray[_moveStack[i]]->getNetList()[j]]->incPartCount(_cellArray[_moveStack[i]]->getPart());
							_netArray[_cellArray[_moveStack[i]]->getNetList()[j]]->decPartCount(!_cellArray[_moveStack[i]]->getPart());
						}
						_partSize[!_cellArray[_moveStack[i]]->getPart()]--;
						_partSize[_cellArray[_moveStack[i]]->getPart()]++;
					}
					for (int i = 0; i <= _bestMoveNum; i++) {
						_cellArray[_moveStack[i]]->setPart(!_cellArray[_moveStack[i]]->getPart());
					}
					for (int i = 0; i < _cellNum; i++)
						_cellArray[i]->unlock();
					_cutSize -= _maxAccGain;
					_bestMoveNum = 0;
					_moveNum = 0;
					_accGain = 0;
					//switch node base on bestMove in _moveStack
				}
				_moveStack.clear();
				break;
			}
			_accGain += _cellArray[_maxGainCell->getId()]->getGain();

			//remove node from bucket
			if (_maxGainCell->getPrev() != NULL)
				_maxGainCell->getPrev()->setNext(_maxGainCell->getNext());
			else
				_bList[_cellArray[_maxGainCell->getId()]->getPart()][_cellArray[_maxGainCell->getId()]->getGain()] = _maxGainCell->getNext();

			if (_maxGainCell->getNext() != NULL)
				_maxGainCell->getNext()->setPrev(_maxGainCell->getPrev());

			//lock the switched node

			_cellArray[_maxGainCell->getId()]->lock();
			_partSize[_cellArray[_maxGainCell->getId()]->getPart()]--;
			_partSize[!_cellArray[_maxGainCell->getId()]->getPart()]++;

			update_gain(_cellArray[_maxGainCell->getId()], max_gain_value);

			if (_accGain > _maxAccGain) {
				_maxAccGain = _accGain;
				_bestMoveNum = _moveNum;
			}
			_moveNum++;
			_moveStack.push_back(_maxGainCell->getId());
		}

	}
	/*****************process end*****************/
}
void Partitioner::update_gain(Cell* target, int& max_gain_value) {
	bool part = target->getPart();
	vector <int> node_need_be_bhanged;
	vector <int>ori_gain;
	max_gain_value = _maxPinNum;

	for (int i = 0; i < target->getNetList().size(); i++) {
		if (_netArray[target->getNetList()[i]]->getPartCount(!part) == 0) {//to size=0
			for (int j = 0; j < _netArray[target->getNetList()[i]]->getCellList().size(); j++) {
					if (find(node_need_be_bhanged.begin(), node_need_be_bhanged.end(), _netArray[target->getNetList()[i]]->getCellList()[j]) == node_need_be_bhanged.end() ) {
						node_need_be_bhanged.push_back(_netArray[target->getNetList()[i]]->getCellList()[j]);
						ori_gain.push_back(_cellArray[_netArray[target->getNetList()[i]]->getCellList()[j]]->getGain());
					}
					_cellArray[_netArray[target->getNetList()[i]]->getCellList()[j]]->incGain();
			}
		}
		 if (_netArray[target->getNetList()[i]]->getPartCount(!part) == 1) {//to size=1
			for (int j = 0; j < _netArray[target->getNetList()[i]]->getCellList().size(); j++) {
				if (_cellArray[_netArray[target->getNetList()[i]]->getCellList()[j]]->getPart() == !part ) {
					if (find(node_need_be_bhanged.begin(), node_need_be_bhanged.end(), _netArray[target->getNetList()[i]]->getCellList()[j]) == node_need_be_bhanged.end() ) {
						node_need_be_bhanged.push_back(_netArray[target->getNetList()[i]]->getCellList()[j]);
						ori_gain.push_back(_cellArray[_netArray[target->getNetList()[i]]->getCellList()[j]]->getGain());
					}
					_cellArray[_netArray[target->getNetList()[i]]->getCellList()[j]]->decGain();
				}
			}
		}
		_netArray[target->getNetList()[i]]->decPartCount(part);
		_netArray[target->getNetList()[i]]->incPartCount(!part);
		if (_netArray[target->getNetList()[i]]->getPartCount(part) == 0) {//from size=0
			for (int j = 0; j < _netArray[target->getNetList()[i]]->getCellList().size(); j++) {
					if (find(node_need_be_bhanged.begin(), node_need_be_bhanged.end(), _netArray[target->getNetList()[i]]->getCellList()[j]) == node_need_be_bhanged.end() ) {
						node_need_be_bhanged.push_back(_netArray[target->getNetList()[i]]->getCellList()[j]);
						ori_gain.push_back(_cellArray[_netArray[target->getNetList()[i]]->getCellList()[j]]->getGain());
					}
					_cellArray[_netArray[target->getNetList()[i]]->getCellList()[j]]->decGain();
			}
		}
		if (_netArray[target->getNetList()[i]]->getPartCount(part) == 1) {//from size=1
			for (int j = 0; j < _netArray[target->getNetList()[i]]->getCellList().size(); j++) {
				if (_cellArray[_netArray[target->getNetList()[i]]->getCellList()[j]]->getPart() == part) {
					if (find(node_need_be_bhanged.begin(), node_need_be_bhanged.end(), _netArray[target->getNetList()[i]]->getCellList()[j]) == node_need_be_bhanged.end() ) {
						node_need_be_bhanged.push_back(_netArray[target->getNetList()[i]]->getCellList()[j]);
						ori_gain.push_back(_cellArray[_netArray[target->getNetList()[i]]->getCellList()[j]]->getGain());
					}
					_cellArray[_netArray[target->getNetList()[i]]->getCellList()[j]]->incGain();
				}
			}
		}
	}
	for (int i = 0; i < node_need_be_bhanged.size(); i++) {
		if (ori_gain[i] != _cellArray[node_need_be_bhanged[i]]->getGain() && !_cellArray[node_need_be_bhanged[i]]->getLock()) {
			if(_cellArray[node_need_be_bhanged[i]]->getNode()->getNext()!=NULL)
				_cellArray[node_need_be_bhanged[i]]->getNode()->getNext()->setPrev(_cellArray[node_need_be_bhanged[i]]->getNode()->getPrev());
			if (_cellArray[node_need_be_bhanged[i]]->getNode()->getPrev() != NULL)
				_cellArray[node_need_be_bhanged[i]]->getNode()->getPrev()->setNext(_cellArray[node_need_be_bhanged[i]]->getNode()->getNext());
			if(_cellArray[node_need_be_bhanged[i]]->getNode()->getPrev()==NULL)
				_bList[_cellArray[node_need_be_bhanged[i]]->getPart()][ori_gain[i]] = _cellArray[node_need_be_bhanged[i]]->getNode()->getNext();
			if(_cellArray[node_need_be_bhanged[i]]->getNode()->getPrev() == NULL && _cellArray[node_need_be_bhanged[i]]->getNode()->getNext() == NULL)
				_bList[_cellArray[node_need_be_bhanged[i]]->getPart()][ori_gain[i]]=NULL;
			_cellArray[node_need_be_bhanged[i]]->getNode()->setNext(NULL);
			_cellArray[node_need_be_bhanged[i]]->getNode()->setPrev(NULL);
			if (_bList[_cellArray[node_need_be_bhanged[i]]->getPart()][_cellArray[node_need_be_bhanged[i]]->getGain()] != NULL) {
				_bList[_cellArray[node_need_be_bhanged[i]]->getPart()][_cellArray[node_need_be_bhanged[i]]->getGain()]->setPrev(_cellArray[node_need_be_bhanged[i]]->getNode());
				_cellArray[node_need_be_bhanged[i]]->getNode()->setNext(_bList[_cellArray[node_need_be_bhanged[i]]->getPart()][_cellArray[node_need_be_bhanged[i]]->getGain()]);
			}
			_bList[_cellArray[node_need_be_bhanged[i]]->getPart()][_cellArray[node_need_be_bhanged[i]]->getGain()] = _cellArray[node_need_be_bhanged[i]]->getNode();
		}
	}
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