//#if defined(__GNUC__)
//#pragma GCC diagnostic push
//#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
//#elif defined(_MSC_VER)
//#pragma warning(disable : 4996)
//#endif
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include<algorithm>
#include<math.h>
#include<set>
#include <ctime>
#include "json.h"

using namespace std;

unordered_map<int,int> idindex;

class relation;
//get a class like above, id and three int
class Inst {
public:
    int id;
    int left;   //number of input ports
    int right;  //number of output ports
    int dual;   //number of dual ports
    int height;
    int r;//行号
    int c;//列号
    pair<int, int> lUp;    //左上角坐标,first is row(x), second is column(y)
    pair<int, int> rDown;  //右下角坐标
    vector<relation> lconnect;
    vector<relation> rconnect;

};

class Stream{
public:
    int id;
    pair<int, int> lidport;
    pair<int, int> ridport;
    vector<pair<int, int> > spath;//actual position of points
};

class relation{
    public:
    relation() {

    }
    Inst* device;
    Stream* stream;
    int port;
};

struct idPort
{
    int id;
    int port;
    bool inOrOut;
};

struct mima {
    int min;
    int max;
};

struct idPortWithStream{
    idPortWithStream()=default;
    idPortWithStream(const pair<int, int>& pair1, const pair<int, int>& pair2, const pair<int, int>& pair3) {
        pointOnBroadWiseLine=pair1;
        pointOnTurn=pair2;
        pointOnDevice=pair3;
    }
    idPortWithStream(const idPortWithStream& idPortWithStream1) {
        pointOnBroadWiseLine=idPortWithStream1.pointOnBroadWiseLine;
        pointOnTurn=idPortWithStream1.pointOnTurn;
        pointOnDevice=idPortWithStream1.pointOnDevice;
    }
    //override operator<< for output
    friend ostream& operator<<(ostream& os, const idPortWithStream& idPortWithStream1) {
        os << "(" << idPortWithStream1.pointOnBroadWiseLine.first << "," << idPortWithStream1.pointOnBroadWiseLine.second << ")";
        os << "(" << idPortWithStream1.pointOnTurn.first << "," << idPortWithStream1.pointOnTurn.second << ")";
        os << "(" << idPortWithStream1.pointOnDevice.first << "," << idPortWithStream1.pointOnDevice.second << ")";
        return os;
    }

    pair<int, int> pointOnBroadWiseLine;
    pair<int, int> pointOnTurn;
    pair<int, int> pointOnDevice;
};

vector<Inst> insts;
vector<Stream> streams;

//read json file to string and ignore blank characters
string readJson(string filename) {
    ifstream in(filename);
    string json;
    string line;
    while (getline(in, line)) {
        json += line;
    }
    return json;
}


//parse the json string like above to class Inst above use jsoncpp
vector<Inst> parseInst(string json) {
    vector<Inst> v;
    Json::Reader reader;
    Json::Value root;

    if (reader.parse(json, root)) {
        for (auto it = root.begin(); it != root.end(); ++it) {
            Inst inst;
            inst.id = stoi(it.key().asString());

            pair<int,int> tmpair(inst.id, distance(root.begin(), it));
            idindex.insert(tmpair);
            inst.left = root[it.key().asString()][0].asInt();
            inst.right = root[it.key().asString()][1].asInt();
            inst.dual = root[it.key().asString()][2].asInt();
            v.push_back(inst);
        }
    }
    return v;
}


//parse the json string like above to class Stream above use jsoncpp
vector<Stream> parseStream(string json) {
    vector<Stream> v;
    Json::Reader reader;
    Json::Value root;
    if (reader.parse(json, root)) {
        for (auto it = root.begin(); it != root.end(); ++it) {
            Stream stream;
            stream.id = it.key().asInt();
            stream.lidport.first = (*it)[0].asInt();
            stream.lidport.second = (*it)[1].asInt();
            stream.ridport.first = (*it)[2].asInt();
            stream.ridport.second = (*it)[3].asInt();
            v.push_back(stream);
        }
    }
    return v;
}

//output the vector<Inst> to a json file using jsoncpp
void outputInst(vector<Inst> v, string filename) {
    Json::Value root;
    Json::FastWriter writer;
    for (auto it = v.begin(); it != v.end(); ++it) {
        root[to_string(it->id)] = Json::Value(Json::arrayValue);
        root[to_string(it->id)][0] = it->lUp.first;
        root[to_string(it->id)][1] = it->lUp.second;
        //root[to_string(it->id)][2] = it->dual;
    }
    ofstream out(filename);
    out << writer.write(root);
}


//find an Inst in the vector of Inst with its id and return the pointer to the Inst
Inst* findInst(vector<Inst>& v, int id) {
//    for (auto it = v.begin(); it != v.end(); ++it) {
//        if (it->id == id) {
//            return &(*it);
//        }
//    }
//    cout<<"NULL"<<endl;
    return &v[idindex[id]];
}

//fill out the lconnect and rconnect of objects int vector<Inst> using vector<Stream>
//use the pair lidport to fill out the rconnect and use the pair reimport to fill out the rconnect
void fillConnect(vector<Inst>& v, vector<Stream>& s) {
    int tmpx=0,tmpy=0;
    for (auto x= v.begin(); x != v.end(); ++x) {
        for (auto& y : s) {
            for (int i = 1; i <= x->left ; ++i) {
                if (x->id==y.ridport.first && i == y.ridport.second) {
                    relation rel;
                    rel.device =findInst(v,y.lidport.first);//y.lidport.first the id of the inst to find
                    rel.stream = &y;
                    rel.port = y.lidport.second;
                    x->lconnect.push_back(rel);
                }
            }
            for (int i = 1; i <= x->right ; ++i) {
                if (x->id==y.lidport.first && i == y.lidport.second) {
                    relation rel;
                    rel.device = findInst(v,y.ridport.first);
                    rel.stream = &y;
                    rel.port = y.ridport.second;
                    x->rconnect.push_back(rel);
                }
            }

        }
        tmpx++;
    }
}

string idPortToString(pair<int,int> idPort,bool inOrOut)
{//make the no duplicate keys
    string currIdPort;
    if(inOrOut)
        currIdPort="In:";
    else
        currIdPort="Out:";
    currIdPort+=to_string(idPort.first);
    currIdPort.push_back(',');
    currIdPort+=to_string(idPort.second);
    return currIdPort;
}

unordered_map<string,int> idPortToStream;//String里面的id是Inst的输入id，int里面的是Stream在vector<Stream>的index
//strings like"Out:33,2"; for no repeated keys; int is the index of the stream in the vector<Stream>
vector<vector<idPort>> streamToInstsPort;//外层vector下标是Stream在vector<Stream>的index，内层vector里面的idPort中的id，是Inst在vector<Inst>的index
//for rewiring outside vector follow the index of the stream maybe lots of starts&ends in one vector<idPort>
vector<vector<int>> streamToInstsId;//外层vector下标是Stream在vector<Stream>的index，内层vector里面是Inst在vector<Inst>的index
//for relocation
vector<vector<int>> instToStreams;//外层vector下标是Inst在vector<Inst>的index，内层vector里面是Stream在vector<Stream>的index
vector<vector<int>> instToInst;//外层vector下标是Inst在vector<Inst>的index，内层vector里面是Inst在vector<Inst>的index


void fillStream(vector<Inst>& v,vector<Stream>& s,int instNum)
{   //all five structures the id means the index

    instToStreams.resize(instNum);
    instToInst.resize(instNum);
    for(auto &stream:s)
    {
        string currStartIDPort=idPortToString(stream.lidport,true);
        string currEndIDPort=idPortToString(stream.ridport,false);
        idPort start;
        start.id=idindex[stream.lidport.first];
        start.port=stream.lidport.second;
        start.inOrOut=true;
        idPort end;
        end.id=idindex[stream.ridport.first];
        end.port=stream.ridport.second;
        end.inOrOut=false;

        int streamId;
        if(idPortToStream.count(currStartIDPort))//make sure each steam has all insts connect to it
            //count() returns the number with a certain key 0/1
        {//if the key is already in the map
            streamId=idPortToStream[currStartIDPort];
            idPortToStream[currEndIDPort]=streamId;             //if there's no key [currEndIDPort] it will creat one
            streamToInstsPort[streamId].push_back(end);
        }
        else if(idPortToStream.count(currEndIDPort))
        {
            streamId=idPortToStream[currEndIDPort];
            idPortToStream[currStartIDPort]=streamId;
            streamToInstsPort[streamId].push_back(start);
        }
        else {
            vector<idPort> tempStreamToInsts;
            tempStreamToInsts.push_back(start);
            tempStreamToInsts.push_back(end);
            streamToInstsPort.push_back(tempStreamToInsts);
            streamId = streamToInstsPort.size() - 1;
            idPortToStream[currStartIDPort] = streamId; //add the start/end string and stream index to the map
            idPortToStream[currEndIDPort] = streamId;
        }
    }
    int streamNum=streamToInstsPort.size();
    streamToInstsId.resize(streamNum);
    for(int i=0;i<streamNum;i++)
    //store id-of-inst vector for each stream (not index of inst)
    {
        //set<int> instSet;
        //vector<int> instSet;
        for(auto &idport:streamToInstsPort[i])
        {
            //instSet.push_back(idindex[idport.id]);
            streamToInstsId[i].push_back(idindex[idport.id]);
        }
//        vector<int> insts;
//        insts.assign(instSet.begin(),instSet.end());//assign():the function fills a %vector with copies of the elements in the range [__first,__last)
//        streamToInstsId.push_back(insts);
    }
    for(int i=0;i<streamNum;i++)
    //push_back the index of insts in each stream in instToStreams
    {
        for(auto &id:streamToInstsId[i])
        {
            instToStreams[id].push_back(i);//id is the index of inst in insts; i is the index of the stream in input
        }
    }

    for(int i=0;i<instNum;i++)
    {
        set<int> instSet;
        for(int j=0;j<instToStreams[i].size();j++)
        {
            int streamIndex=instToStreams[i][j];//to every stream
            for(int k=0;k<streamToInstsId[streamIndex].size();k++)
            {
                instSet.insert(streamToInstsId[streamIndex][k]);//include the inst itself, what about the loop?
            }
        }
        instToInst[i].assign(instSet.begin(),instSet.end());
    }

//    for(auto it=idPortToStream.begin(); it != idPortToStream.end(); ++it){
//        cout<<it->first<<" "<<it->second<<endl;
//    }
}

//模拟退火调整列间布局
//rCnt是元件排列后行数
//cCnt是元件排列后列数
//int** inst2d是记录i行j列元件编号的二维数组(instance two dimension)
void SA_relocation(vector<Inst> insts, vector<Stream> streams,int cCnt, vector<int> vecArr[])
{

    //设置模拟退火参数
    const double max_temper = 1000, min_temper = 1;   //最高最低温度 1000
    int iterT = 1000;  //迭代次数 1000
    double temprature = max_temper;
    double dec = 0.99;

    //记录模拟退火的优化程度
    double grow = 0;
    double temgrow = 0;

    //优化目标为线网的跨度较小，即跨越的行较少为目标
    //每次随机选取一列，并在列内随机选择两个元素进行交换
    srand((int)time(0));
    int cRan = 0;  //column random
    int rRan1 = 0, rRan2 = 0;
    while (temprature > min_temper)
    {
        for (int i = 0; i < iterT; i++)
        {
            temgrow = 0;
            cRan = rand() % cCnt;
            if (vecArr[cRan].size() == 1)
                continue;
            rRan1 = rand() % vecArr[cRan].size();
            rRan2 = rand() % vecArr[cRan].size();
            while (rRan1 == rRan2) {
                rRan1 = rand() % vecArr[cRan].size();
                rRan2 = rand() % vecArr[cRan].size();
            }

            //计算在交换两个元件位置之后，线网的宽度改变

            //考虑与第一个元件连接的线，计算它们的线网宽度改变
            //左边
            for (int i = 0; i < insts[vecArr[cRan][rRan1]].left; i++)//i is port number on the left
            {
                if(vecArr[cRan][rRan1]==-1)
                    break;//empty inst, only consider the other one
                pair<int, int> tmppair(insts[vecArr[cRan][rRan1]].id,i);
                auto pair=idPortToStream.find(idPortToString(tmppair,true));
                if(pair==idPortToStream.end()) continue;
                int stream = pair->second;
                int old = abs(insts[streams[stream].lidport.first].r - rRan1);
                int now = abs(insts[streams[stream].lidport.first].r - rRan2);
                temgrow += now - old;
            }
            //右边
            for (int i = 0; i < insts[vecArr[cRan][rRan1]].right; i++)
            {
                if(vecArr[cRan][rRan1]==-1)
                    break;//empty inst, only consider the other one
                pair<int, int> tmppair(insts[vecArr[cRan][rRan1]].id,i);
                auto pair=idPortToStream.find(idPortToString(tmppair,false));
                if(pair==idPortToStream.end()) continue;
                int stream = pair->second;
                int old = abs(insts[streams[stream].ridport.first].r - rRan1);
                int now = abs(insts[streams[stream].ridport.first].r - rRan2);
                temgrow += now - old;
            }
            //考虑与第二个元件连接的线，计算它们的线网宽度改变
            //左边
            for (int i = 0; i < insts[vecArr[cRan][rRan2]].left; i++)
            {
                if(vecArr[cRan][rRan2]==-1)
                    break;//empty inst, only consider the other one
                pair<int, int> tmppair(insts[vecArr[cRan][rRan2]].id,i);
                auto pair=idPortToStream.find(idPortToString(tmppair,true));
                if(pair==idPortToStream.end()) continue;
                int stream = pair->second;
                int old = abs(insts[streams[stream].lidport.first].r - rRan2);
                int now = abs(insts[streams[stream].lidport.first].r - rRan1);
                temgrow += now - old;
            }
            //右边
            for (int i = 0; i < insts[vecArr[cRan][rRan2]].right; i++)
            {
                if(vecArr[cRan][rRan2]==-1)
                    break;//empty inst, only consider the other one
                pair<int, int> tmppair(insts[vecArr[cRan][rRan2]].id,i);
                auto pair=idPortToStream.find(idPortToString(tmppair,false));
                if(pair==idPortToStream.end()) continue;
                int stream = pair->second;
                int old = abs(insts[streams[stream].ridport.first].r - rRan2);
                int now = abs(insts[streams[stream].ridport.first].r - rRan1);
                temgrow += now - old;
            }
            //判断是否需要接受当前的交换操作
            if (temgrow < 0 || (rand() % 1000) < exp((-temgrow) / temprature) * 1000)
            {
                int temp = vecArr[cRan][rRan1];
                vecArr[cRan][rRan1] = vecArr[cRan][rRan2];
                vecArr[cRan][rRan2] = temp;
                grow += temgrow;
                temprature = temprature * dec;
                //cout<<"grow:"<<grow<<"  "<<"temp: "<<temprature<<endl;
            }
        }
        //temprature = temprature * dec;
    }
}

//compare the inst.left
bool cmp1(Inst a,Inst b){
    return a.left>b.left;
}
//compare the inst.right
bool cmp2(Inst a,Inst b){
    return a.right<b.right;
}

vector<int> cntConnect;

bool cmp3(Inst a,Inst b)//,vector<int> cnt)
{
    if (cntConnect[idindex.find(a.id)->second]<cntConnect[idindex.find(b.id)->second])
        return true;
    else
        return false;
}

void countConnection(const vector<int> &inport, vector<int> &cntConnect) {
    for (int i = 0; i < inport.size(); i++) {
        //for (int j = 0; j < size; j++) {
            for(int k=0;k<instToInst[inport[i]].size();k++) {
                if(::cntConnect[instToInst[inport[i]][k]] == -1) continue;
                ::cntConnect[instToInst[inport[i]][k]]++;
            }
        //}

    }
}

vector<Inst> &
getNextCol(const vector<Inst> &insts, const vector<int> &inport, const vector<int> &outport,
           vector<Inst> &instsTemp, vector<int> &cntConnect,int& Major){

    if(inport.size() < outport.size()){
        Major=0;
        //size=inport.size();
        countConnection(inport, cntConnect);
    }else{
        Major=1;
        //size=outport.size();
        countConnection(outport, cntConnect);
    }
    for(int i = 0; i < cntConnect.size(); i++){
        if(cntConnect[i]>0){
            instsTemp.push_back(insts[i]);
        }
    }
    sort(instsTemp.begin(),instsTemp.end(),cmp3);
    return ::insts;
}

void relocation(vector<Inst>& insts, vector<Stream>& streams, vector<int> vecArr[], vector<int>& widthVec,
                int& vecNum, int& totalWidth, int& tmpCol, int& maxLine) {
    //inst.height=max(inst.right,inst.left)+1, compute height for each inst
    for (auto &inst: insts) {
        inst.height = 2 * max(inst.right, inst.left) + 1;
    }

    vector<int> inport;//find every in/out port
    vector<int> outport;
    cntConnect.resize(insts.size());
    for (int i = 0; i < insts.size(); i++) {
        if (insts[i].left == 0) {
            inport.push_back(i);
            cntConnect[i] = -1;
        } else if (insts[i].right == 0) {
            outport.push_back(i);
            cntConnect[i] = -1;
        } else {
            cntConnect[i] = 0;
        }
    }


    //计算左端口多还是右端口多
    //int Major = 0;//0:left,1:right
    //int size;

    vector<Inst> instsTemp;//store the insts after comparison

    /*
    if (inport.size() > outport.size()) {
        Major = 0;
        size = inport.size();
        countConnection(inport, cntConnect);
    } else {
        Major = 1;
        size = outport.size();
        countConnection(outport, cntConnect);
    }
    */
    //as default: follow inport
    countConnection(inport, cntConnect);
    for (int i = 0; i < cntConnect.size(); i++) {
        if (cntConnect[i] > 0) {
            instsTemp.push_back(insts[i]);
        }
    }

    sort(instsTemp.begin(), instsTemp.end(), cmp3);
    /*
    int lc = 0, rc = 0;
    for (int i = 0; i < insts.size(); i++)
    {
        lc += insts[i].left;
        rc += insts[i].right;
    }
    vector<Inst> instsTemp;//store the insts after comparison

    if (lc >= rc)
    {//以左边为基准
        instsTemp= insts;
        sort(instsTemp.begin(),instsTemp.end(),cmp1);
        Major=0;
    }
    else
    {//以右边为基准
        instsTemp = insts;
        sort(instsTemp.begin(), instsTemp.end(), cmp2);
        Major=1;
    }
    */

    //compute the max height of insts
    int maxInstHeight = 0, connectEndPoint = 0;
    for (int i = 0; i < insts.size(); i++) {
        maxInstHeight = maxInstHeight > insts[i].height ? maxInstHeight : insts[i].height;
        //if(cntConnect[idindex[instsTemp[i].id]]>0) connectEndPoint=i;
    }
    //max是一列元件占用空间的上限，maxNum是列内元件数的上限
    maxInstHeight = maxInstHeight * ceil(sqrt(double(insts.size())));
    int maxNum = ceil(sqrt(double(insts.size())));

    //pre-relocation
    //p存放已经确定位置的全部元件序号
    //vector<int> p;
    int p=0;

    //记录每一列占据的纵向网格数，colVec[i]，表示第i列的纵向宽度
    vector<int> colVec;
    //记录每一列外层包围的空间厚度，maxCrowdVec[i]表示第i列的横向宽度
    vector<int> maxCrowdVec;

    //记录当前元件和前面所有列内元件连线数的总和
    int sum = 0;//total stream number of insts in the column
    int instNum = 0;
    int colMax = 0;//save the maxium height of columns


    for (int i = 0; i < maxNum; i++) {
        if (i < inport.size()) {
            vecArr[0].push_back(inport[i]);
            //else
            //vecArr[0].push_back(-1);
            //p.push_back(insts[inport[i]].id);
            p++;
        }
    }
    vecNum++;
    while (true) {
        //compute the total stream&insts number in the column
        sum = 0;
        for (instNum = 0;
             sum < maxInstHeight    //space limitation
             && instNum < instsTemp.size()
             && instNum < maxNum          // no too much insts in one stream
            //&& instNum<connectEndPoint
                ;
             instNum++) {
            sum += instToStreams[idindex[instsTemp[instNum].id]].size();
        }
        if (sum > maxInstHeight) {//最后一个加多了再吐一个，上面的循环因为sum>=max终止时就需要这样做
            sum -= instToStreams[idindex[instsTemp[instNum].id]].size();
            instNum--;
        }
        //if(instNum >= connectEndPoint) connectEndPoint==INT_MAX;
        int temWidth = 0, temBigWidth = 0, thick = 0, temInd = 0, temCol = 0;//, streamsPreCol = 0;
        //add number of instNum instances to the column vecArr and compute the width and height of the column
        for (int i = 0; i < instNum; i++) {
            int id = instsTemp[i].id;
            vecArr[vecNum].push_back(idindex[id]);//加入靠前的device,记录index
            cntConnect[idindex[id]]=-1;
            //instsTemp.erase(instsTemp.begin() + i);
            thick += instToStreams[id].size();//根据所连的线的总数，更新需要和下一列留出多少间隔

            //p.push_back(id);
            p++;
            if (temBigWidth < instToStreams[id].size()) {
                temBigWidth = instToStreams[id].size();//更新单个元件最大连线数
            }
            temCol += (3 * instToStreams[id].size() + 1);//更新这一列全局高度,3*连线数+1是为了预留空间，原来是两倍连线数
        }

        widthVec.push_back(thick * 2 + 4 + temBigWidth);//记录每一列占据的横向网格数
        colVec.push_back(temCol);
        if (temCol > colMax)
            colMax = temCol;
        totalWidth += (thick * 2 + 4 + temBigWidth);
        maxCrowdVec.push_back(thick + 2);
        //vecNum++;//总列数

        if (p+outport.size() == insts.size()) {
            instsTemp.clear();
            vector<Inst>().swap(instsTemp);//与空临时变量交换，释放内存
            vecNum++;
            break;//已经全部放入，此为唯一的退出条件
        } else {
//                auto ins=instsTemp.begin();
//                instsTemp.erase(ins,ins+instNum);
//                connectEndPoint-=instNum;
            //cntConnect.clear();
            countConnection(vecArr[vecNum++], cntConnect);
            instsTemp.clear();
            vector<Inst>().swap(instsTemp);
            for (int i = 0; i < cntConnect.size(); i++) {
                if (cntConnect[i] > 0) {
                    instsTemp.push_back(insts[i]);
                }
            }
            sort(instsTemp.begin(), instsTemp.end(), cmp3);
        }

    }

    for (int i = 0; i < maxNum; i++) {
        if (i < outport.size()) {
            vecArr[vecNum].push_back(outport[i]);
            //p.push_back(insts[outport[i]].id);
        }
    }
    vecNum++;
    //end of the prerelocation

    //按元件序号的顺序记录每个元件的行列位置
    for (int i = 0; i < vecNum; i++) {//i cols，j rows
        for (int j = 0; j < vecArr[i].size(); j++) {
            insts[vecArr[i][j]].r=j;
            insts[vecArr[i][j]].c=i;
        }
    }

    //模拟退火
    SA_relocation(insts,streams,vecNum,vecArr);
    //退火后，按元件序号的顺序重新记录每个元件的行列位置
    for (int i = 0; i < vecNum; i++) {//i cols，j rows
        for (int j = 0; j < vecArr[i].size(); j++) {
            insts[vecArr[i][j]].r=j;
            insts[vecArr[i][j]].c=i;
        }
        if(vecArr[i].size()>maxLine) maxLine=vecArr[i].size();
    }


    //确定元件的物理坐标
    int pin = 2, idWidth = 8;//the width of the insts
    int rowMaxNum = 0;   //最大行数
    for (int i = 0; i < vecNum; i++)
    {
        if (vecArr[i].size() > rowMaxNum) rowMaxNum = vecArr[i].size();
    }
    int* width = new int[vecNum];   //存储每列元件左上角位置横向值
    int cnt = 0;//2;
    for (int i = 0; i < vecNum; i++)  //列
    {
        width[i] = cnt;
        for (int j = 0; j < vecArr[i].size(); j++)  //此列右端out端口数
        {
            cnt += insts[vecArr[i][j]].right;
        }
        for (int j = 0; i != vecNum - 1 && j < vecArr[i + 1].size(); j++)  //下一列左端in端口数，若为最后一列则不考虑
        {
            cnt += insts[vecArr[i+1][j]].left;
        }
        cnt += pin * 2 + idWidth;//two cols of pins and one col of insts
    }
    //计算每个元件左上角
    int h = 0;
    for (int i = 0; i < rowMaxNum; i++)  //行
    {
        int hGap = streamToInstsPort.size() / vecArr[i].size();
        int maxHeight = 0;
        for (int j = 0; j < vecNum; j++)   //列
        {
            if (i >= vecArr[j].size() || vecArr[j][i] == -1)  //此行没有此列元件,跳过
                continue;
            insts[vecArr[j][i]].height = 2 * (max(insts[vecArr[j][i]].left, insts[vecArr[j][i]].right) + insts[vecArr[j][i]].dual) + 1;
            insts[vecArr[j][i]].lUp.first = width[j];
            insts[vecArr[j][i]].lUp.second = h;
            insts[vecArr[j][i]].rDown.first = width[j] + idWidth;
            insts[vecArr[j][i]].rDown.second = h + insts[vecArr[j][i]].height;
            if (insts[vecArr[j][i]].height > maxHeight) maxHeight = insts[vecArr[j][i]].height;
        }
        h += hGap + maxHeight;
    }


}


unordered_map<string,idPortWithStream> idPortWithStreamRes;

/*
void outputStream(string filename) {
    Json::Value root;
    Json::FastWriter writer;
    for (auto stream:streams) {
        //root[to_string(stream.id)] = Json::Value(Json::arrayValue);
        int inputStartId=stream.lidport.first;
        int inputEndId=stream.ridport.first;
        int actualStartId=idindex[inputStartId];
        int actualEndId=idindex[inputEndId];
        pair<int,int> startIdPort(actualStartId,stream.lidport.second);
        pair<int,int> endIdPort(actualEndId,stream.ridport.second);
        //pair<int,int> startIdPort(stream.lidport.first,stream.lidport.second);
        //pair<int,int> endIdPort(stream.ridport.first,stream.ridport.second);
        string startIdPortStr=idPortToString(startIdPort,true);
        string endIdPortStr=idPortToString(endIdPort,false);
        idPortWithStream startIdPortWithStream=idPortWithStreamRes[startIdPortStr];
        idPortWithStream endIdPortWithStream=idPortWithStreamRes[endIdPortStr];
        string tmp;
        //name
        tmp.append(to_string(stream.lidport.first));
        tmp.append(" ");
        tmp.append(to_string(stream.lidport.second));
        tmp.append(" ");
        tmp.append(to_string(stream.ridport.first));
        tmp.append(" ");
        tmp.append(to_string(stream.ridport.second));
        string name=tmp;
        root[name]  = Json::Value(Json::arrayValue);
        tmp.clear();
        //line board
        tmp.append(to_string(startIdPortWithStream.pointOnBroadWiseLine.first));
        tmp.append(" ");
        tmp.append(to_string(startIdPortWithStream.pointOnBroadWiseLine.second));
        tmp.append(" ");
        tmp.append(to_string(endIdPortWithStream.pointOnBroadWiseLine.first));
        tmp.append(" ");
        tmp.append(to_string(endIdPortWithStream.pointOnBroadWiseLine.second));
        root[name][0] = tmp;
        tmp.clear();
        //line turn for start
        tmp.append(to_string(startIdPortWithStream.pointOnBroadWiseLine.first));
        tmp.append(" ");
        tmp.append(to_string(startIdPortWithStream.pointOnBroadWiseLine.second));
        tmp.append(" ");
        tmp.append(to_string(startIdPortWithStream.pointOnTurn.first));
        tmp.append(" ");
        tmp.append(to_string(startIdPortWithStream.pointOnTurn.second));
        root[name][1] = tmp;
        tmp.clear();
        //line turn for end
        tmp.append(to_string(endIdPortWithStream.pointOnBroadWiseLine.first));
        tmp.append(" ");
        tmp.append(to_string(endIdPortWithStream.pointOnBroadWiseLine.second));
        tmp.append(" ");
        tmp.append(to_string(endIdPortWithStream.pointOnTurn.first));
        tmp.append(" ");
        tmp.append(to_string(endIdPortWithStream.pointOnTurn.second));
        root[name][2] = tmp;
        tmp.clear();

        //line device for start
        tmp.append(to_string(startIdPortWithStream.pointOnTurn.first));
        tmp.append(" ");
        tmp.append(to_string(startIdPortWithStream.pointOnTurn.second));
        tmp.append(" ");
        tmp.append(to_string(startIdPortWithStream.pointOnDevice.first));
        tmp.append(" ");
        tmp.append(to_string(startIdPortWithStream.pointOnDevice.second));
        root[name][3] = tmp;
        tmp.clear();
        //line turn for end
        tmp.append(to_string(endIdPortWithStream.pointOnTurn.first));
        tmp.append(" ");
        tmp.append(to_string(endIdPortWithStream.pointOnTurn.second));
        tmp.append(" ");
        tmp.append(to_string(endIdPortWithStream.pointOnDevice.first));
        tmp.append(" ");
        tmp.append(to_string(endIdPortWithStream.pointOnDevice.second));
        root[name][4] = tmp;
        tmp.clear();
    }
    ofstream out(filename);
    out << writer.write(root);
}
*/
void outputStream(string filename) {
    Json::Value root;
    Json::FastWriter writer;
    for (auto stream:streams) {
        //root[to_string(stream.id)] = Json::Value(Json::arrayValue);
        int inputStartId=stream.lidport.first;
        int inputEndId=stream.ridport.first;
        int actualStartId=idindex[inputStartId];
        int actualEndId=idindex[inputEndId];
        pair<int,int> startIdPort(actualStartId,stream.lidport.second);
        pair<int,int> endIdPort(actualEndId,stream.ridport.second);
        //pair<int,int> startIdPort(stream.lidport.first,stream.lidport.second);
        //pair<int,int> endIdPort(stream.ridport.first,stream.ridport.second);
        string startIdPortStr=idPortToString(startIdPort,true);
        string endIdPortStr=idPortToString(endIdPort,false);
        idPortWithStream startIdPortWithStream=idPortWithStreamRes[startIdPortStr];
        idPortWithStream endIdPortWithStream=idPortWithStreamRes[endIdPortStr];
        string tmp;
        //name
        tmp.append(to_string(stream.lidport.first));
        tmp.append(" ");
        tmp.append(to_string(stream.lidport.second));
        tmp.append(" ");
        tmp.append(to_string(stream.ridport.first));
        tmp.append(" ");
        tmp.append(to_string(stream.ridport.second));
        string name=tmp;
        root[name]  = Json::Value(Json::arrayValue);
        tmp.clear();
        //line device for start
        tmp.append(to_string(startIdPortWithStream.pointOnDevice.first));
        tmp.append(" ");
        tmp.append(to_string(startIdPortWithStream.pointOnDevice.second));
        tmp.append(" ");
        tmp.append(to_string(startIdPortWithStream.pointOnTurn.first));
        tmp.append(" ");
        tmp.append(to_string(startIdPortWithStream.pointOnTurn.second));
        root[name][0] = tmp;
        tmp.clear();
        //line turn for start
        tmp.append(to_string(startIdPortWithStream.pointOnTurn.first));
        tmp.append(" ");
        tmp.append(to_string(startIdPortWithStream.pointOnTurn.second));
        tmp.append(" ");
        tmp.append(to_string(startIdPortWithStream.pointOnBroadWiseLine.first));
        tmp.append(" ");
        tmp.append(to_string(startIdPortWithStream.pointOnBroadWiseLine.second));
        root[name][1] = tmp;
        tmp.clear();
        //line board
        tmp.append(to_string(startIdPortWithStream.pointOnBroadWiseLine.first));
        tmp.append(" ");
        tmp.append(to_string(startIdPortWithStream.pointOnBroadWiseLine.second));
        tmp.append(" ");
        tmp.append(to_string(endIdPortWithStream.pointOnBroadWiseLine.first));
        tmp.append(" ");
        tmp.append(to_string(endIdPortWithStream.pointOnBroadWiseLine.second));
        root[name][2] = tmp;
        tmp.clear();

        //line turn for end
        tmp.append(to_string(endIdPortWithStream.pointOnBroadWiseLine.first));
        tmp.append(" ");
        tmp.append(to_string(endIdPortWithStream.pointOnBroadWiseLine.second));
        tmp.append(" ");
        tmp.append(to_string(endIdPortWithStream.pointOnTurn.first));
        tmp.append(" ");
        tmp.append(to_string(endIdPortWithStream.pointOnTurn.second));
        root[name][3] = tmp;
        tmp.clear();
        //line turn for end
        tmp.append(to_string(endIdPortWithStream.pointOnTurn.first));
        tmp.append(" ");
        tmp.append(to_string(endIdPortWithStream.pointOnTurn.second));
        tmp.append(" ");
        tmp.append(to_string(endIdPortWithStream.pointOnDevice.first));
        tmp.append(" ");
        tmp.append(to_string(endIdPortWithStream.pointOnDevice.second));
        root[name][4] = tmp;
        tmp.clear();
    }
    ofstream out(filename);
    out << writer.write(root);
}

void printLines()//unordered_map<string,idPortWithStream> idPortWithStreamRes,vector<Stream> streams,unordered_map<int,int> idindex)
{
    for(auto stream:streams)
    {
        int inputStartId=stream.lidport.first;
        int inputEndId=stream.ridport.first;
        int actualStartId=idindex[inputStartId];
        int actualEndId=idindex[inputEndId];
//        pair<int,int> startIdPort(actualStartId,stream.lidport.second);
//        pair<int,int> endIdPort(actualEndId,stream.ridport.second);
        pair<int,int> startIdPort(inputStartId,stream.lidport.second);
        pair<int,int> endIdPort(inputEndId,stream.ridport.second);
        string startIdPortStr=idPortToString(startIdPort,true);
        string endIdPortStr=idPortToString(endIdPort,false);
        idPortWithStream startIdPortWithStream=idPortWithStreamRes[startIdPortStr];
        idPortWithStream endIdPortWithStream=idPortWithStreamRes[endIdPortStr];
        //line1:startIdPortWithStream.pointOnBroadWiseLine,endIdPortWithStream.pointOnBroadWiseLine;


    }
}

typedef pair<int,int> colwid;
// unordered_map<string,int> net::newendwiseWire(int streamId,vector<int> vecArr[], vector<devMatric> dmVec[], int vecNum,int**net,int *xfull,mima gapX[],vector<int> devStrNum[],int *yfull,mima gapY[],bool FirstOrNot,bool LastOrNot,int broadWireX,int currentColumn,vector<int> currentColDev,int inOrOut,int broadWireGap,vector<int> devUpFull[],int streamIndex,bool mossSecondWire)
//inOrout：当前列元件与当前线连接的是左端口还是右端口
//currentColumn：当前列号
//idPortForCurrStream：当前线连接的idPort，仅有左端口或右端口
//idPortWithStreamRes：按元件号端口号inout进行索引，存储当前端口对应的三个点
//broadWireX：当前线的x坐标
//函数得到端口连接到横向线的三个点

void newendwiseWire(bool inOrOut,int currentColumn,vector<idPort> idPortForCurrStream,vector<Inst> insts,
                    //unordered_map<string,idPortWithStream> &idPortWithStreamRes,
                    int broadWireX,mima gapX[],double gapXCrowd[],mima gapY[],int *xfull,int* yfull)
{


    //计算纵向线y坐标及其所在列间隙，in则计算左侧间隙，out计算右侧间隙
    int endWireY=0;
    int endWireGap=-1;
    //记录当前端口的y坐标
    int devY=0;
    if(inOrOut)
    {
        devY=gapY[currentColumn].max-1;//+1;
        endWireY=gapY[currentColumn].max-yfull[currentColumn*2];
        endWireGap=currentColumn*2;
        yfull[endWireGap]++;
    }
    else
    {
        devY=gapY[currentColumn+1].min+1;//-1;
        endWireY=gapY[currentColumn+1].min+yfull[currentColumn*2+1];
        endWireGap=currentColumn*2+1;
        yfull[endWireGap]++;
    }


    //绘制纵向线和横向短线
    //确定当前列与横向线连接点
    pair<int,int> pointOnBroadWiseLine(endWireY,broadWireX);


    for(int i=0;i<idPortForCurrStream.size();i++)
    {
        idPort currIdPort=idPortForCurrStream[i];
        string idPortStr=idPortToString(pair<int,int>(currIdPort.id,currIdPort.port),!inOrOut);
        int devX=insts[currIdPort.id].lUp.second+currIdPort.port*2-1;   //start at 1 distance from top
        pair<int,int> pointOnDevice(devY,devX);
        pair<int,int> pointOnTurn(endWireY,devX);
        idPortWithStream currIdPortWithStream(pointOnBroadWiseLine,pointOnTurn,pointOnDevice);
        idPortWithStreamRes[idPortStr]=currIdPortWithStream;
    }
}




// void  net::newbroadwiseWire(vector<int> vecArr[],vector<devMatric> dmVec[], int vecNum, int maxLine,int**net,int *xfull,mima gapX[],double gapXCrowd[],int *yfull,mima gapY[],vector<colwid> devForStream,vector<int> devStrNum[],int inOrOut,int width,vector<int> devUpFull[],int streamIndex)
//streamIndex：当前线号
//devForStream：当前线连接的元件的行列号，按列号从大到小排序，列号一样按行号从大到小排序
//idPortForCurrStream：当前线连接的idPort，按列排序
void  newbroadwiseWire(vector<int> vecArr[],vector<Inst>& insts,int streamIndex,vector<colwid> devForStream,vector<idPort> idPortForCurrStream,
                       mima gapX[],double gapXCrowd[],mima gapY[],int *xfull,int* yfull,int& maxLine)
{
    //计算和当前线有连接的各列元件的最小和最大行号中距离最小的最小行号和最大行号
    //得到横向线的布线最佳范围
    //需要修改dmVec部分
    int minHighLine=maxLine,maxLowLine=0;
    //对应的行间隙坐标
    int minHighLineX=0,maxLowLineX=0;
    int tempDevIndex=0;
    while(tempDevIndex<devForStream.size())
    {
        int currentColumn=devForStream[tempDevIndex].second;
        if(devForStream[tempDevIndex].first<minHighLine)
        {
            minHighLine=devForStream[tempDevIndex].first;
            minHighLineX=insts[vecArr[currentColumn][minHighLine]].lUp.second;
        }

        int LowLine;
        while(tempDevIndex<devForStream.size()&&currentColumn==devForStream[tempDevIndex].second)
        {
            LowLine=devForStream[tempDevIndex].first;
            tempDevIndex++;
        }
        if(LowLine>maxLowLine)
        {
            maxLowLine=LowLine;
            maxLowLineX=insts[vecArr[currentColumn][maxLowLine]].rDown.second;
        }

    }

    //获取最佳的横向线布线行间隙
    int broadWireX=-1;
    int broadWireGap=-1;
    double minNum=10;
    //如果在最佳范围内，存在没有满的行间隙，则选择拥挤度最低的行间隙
    for(int currentGap=0;currentGap<=maxLine;currentGap++)
    {
        if(gapX[currentGap].max<minHighLineX||gapX[currentGap].min>maxLowLineX||gapXCrowd[currentGap]>=1)
            continue;
        if(gapXCrowd[currentGap]<minNum)
        {
            minNum=gapXCrowd[currentGap];
            broadWireGap=currentGap;
        }
    }
    //如果在最佳范围内，不存在没有满的行间隙，则选择离最佳范围最近的行间隙
    if(minNum==10)
    {
        for(int currentGap=0;currentGap<=maxLine;currentGap++)
        {
            if(gapXCrowd[currentGap]>=1)
                continue;
            double influenceNum=abs(gapX[currentGap].min+xfull[currentGap]-1.0*(minHighLineX+maxLowLineX)/2)/(gapX[currentGap].min+xfull[currentGap]);
            if(minNum>gapXCrowd[currentGap]+influenceNum)
            {
                minNum=gapXCrowd[currentGap]+influenceNum;
                broadWireGap=currentGap;
            }
        }
    }
    //如果最佳范围为一行，即选择该行的上侧未满行间隙或下侧未满行间隙
    if(minHighLine==maxLowLine)
    {
        if(gapXCrowd[minHighLine+1]<1)
            broadWireGap=minHighLine+1;
        else if(gapXCrowd[minHighLine]<1)
            broadWireGap=minHighLine;
    }


    //计算横向线布线位置
    broadWireX=gapX[broadWireGap].min+xfull[broadWireGap];
    if(broadWireGap==0)
        broadWireX=gapX[broadWireGap].max-xfull[broadWireGap];
    xfull[broadWireGap]++;
    gapXCrowd[broadWireGap]=1.0*xfull[broadWireGap]/(gapX[broadWireGap].max-gapX[broadWireGap].min);

    //对每一列元件进行处理，调用纵向布线函数
    int currentDevIndex=0;

    while(currentDevIndex<idPortForCurrStream.size())
    {
        idPort currIdPort=idPortForCurrStream[currentDevIndex];
        int currentColumn=insts[currIdPort.id].c;

        vector<idPort> inIdPort;
        vector<idPort> outIdPort;
        //把每一列元件放至容器中
        while(currentDevIndex<idPortForCurrStream.size())//&&insts[currIdPort.id].c==currentColumn)
        {
            currIdPort=idPortForCurrStream[currentDevIndex];
            if(insts[currIdPort.id].c!=currentColumn) break;
            //if(!idPortForCurrStream[currentDevIndex].inOrOut)
            if(!currIdPort.inOrOut)
                inIdPort.push_back(currIdPort);
            else
                outIdPort.push_back(currIdPort);
            currentDevIndex++;
            //currIdPort=idPortForCurrStream[currentDevIndex];
        }

        if(inIdPort.size()>0)
            newendwiseWire(true,currentColumn,inIdPort,insts,broadWireX,gapX,gapXCrowd,gapY,xfull,yfull);
        if(outIdPort.size()>0)
            newendwiseWire(false,currentColumn,outIdPort,insts,broadWireX,gapX,gapXCrowd,gapY,xfull,yfull);
    }
    //printLines(idPortWithStreamRes,streams,idindex);
}

bool cmpIdPortForStreamByCol(idPort p1,idPort p2)
{
    return insts[idindex[p1.id]].c>insts[idindex[p2.id]].c;
}
bool cmpColWid(colwid cw1,colwid cw2)
{
    if(cw1.second==cw2.second)
        return cw1.first>cw2.first;
    return cw1.second>cw2.second;
}
//布线函数
void  newWiring(vector<int> vecArr[],
                     //vector<devMatric> dmVec[], 
                     int& vecNum, int& temCol, int& width,int& maxLine)
{
    //行间隙中已经被占用的空间
    int xfull[maxLine+1];
    //行间隙的上下界
    mima gapX[maxLine+1];
    //行间隙的拥挤度
    double gapXCrowd[maxLine+1];
    int previousLineLow=0;
    //计算行间隙的上下界，初始化行间隙被占用的空间，及行间隙的拥挤度
    for(int j=0;j<maxLine;j++)
    {
        int nowGapLow=0;
        int nowGapHigh=0;
        nowGapHigh=previousLineLow;
        for(int i=0;i<vecNum;i++)
        {
            if(j>=vecArr[i].size())
                continue;
            nowGapLow=insts[vecArr[i][j]].lUp.second-1;//-1
            //nowGapLow=insts[vecArr[i][j]].rDown.second-1;
            if(insts[vecArr[i][j]].rDown.second+1>previousLineLow)
                previousLineLow=insts[vecArr[i][j]].rDown.second+1;
        }
        mima nowGapLH={nowGapHigh,nowGapLow};
        gapX[j]=nowGapLH;
        if(j==maxLine-1)
            gapX[maxLine].min=previousLineLow;
        gapXCrowd[j]=0;
        xfull[j]=1;
    }
    xfull[maxLine]=1;
    gapX[maxLine].max=temCol-1;


    //yfull[i*2]为第i列的左侧列间隙占用量,yfull[i*2+1]为第i列的右侧列间隙占用量
    int yfull[vecNum*2+2];
    //列间隙的左右界
    mima gapY[vecNum+1];
    int previousColumnRight=0;
    //初始化列间隙占用量和列间隙左右界
    for(int i=0;i<vecNum;i++)
    {
        int nowGapLeft=0;
        int nowGapRight=0;
        nowGapLeft=previousColumnRight;
        for(int j=0;j<vecArr[i].size();j++)
        {
            nowGapRight=insts[vecArr[i][j]].lUp.first-1;
            if(insts[vecArr[i][j]].rDown.first+1>previousColumnRight)
                previousColumnRight=insts[vecArr[i][j]].rDown.first+1;
        }
        mima nowGapLR={nowGapLeft,nowGapRight};
        gapY[i]=nowGapLR;
        if(i==vecNum-1)
            gapY[vecNum].min=previousColumnRight;
        yfull[i*2]=1;
        yfull[i*2+1]=1;
    }
    yfull[vecNum*2]=1;
    yfull[vecNum*2+1]=1;

    //对所有线进行下一步的操作
    for(int i=0;i<streamToInstsPort.size();i++)
    {
        int currentStream=i;
        //记录当前线的连接元件的行列坐标，并排序
        vector<colwid> devForStream;
        for(int j=0;j<streamToInstsId[i].size();j++)
        {
            int instId=insts[streamToInstsId[i][j]].id;
            int col=insts[instId].c;
            int row=insts[instId].r;
            colwid currentDevCR(row,col);
            devForStream.push_back(currentDevCR);
        }
        sort(devForStream.begin(),devForStream.end(),cmpColWid);
        //vector<int> idPortForCurrStream=streamToInstsId[i];
        vector<idPort> idPortForCurrStream=streamToInstsPort[i];    //the idPort in it:id is index
        sort(idPortForCurrStream.begin(),idPortForCurrStream.end(),cmpIdPortForStreamByCol);

        newbroadwiseWire(vecArr,insts,i,devForStream,idPortForCurrStream,
                gapX,gapXCrowd,gapY,xfull,yfull,maxLine);

    }
    return;
}



int main() {
    string json = readJson("../inst.json");
    insts = parseInst(json);

    json=readJson("../net.json");
    streams = parseStream(json);

    fillConnect(insts, streams);

    fillStream(insts, streams,insts.size());

    //output the map idindex
//    for(auto key:idindex){
//        cout<<key.first<<' '<<key.second<<endl;
//    }


    //int Max_cols=sqrt(insts.size());
    vector<int> vecArr[insts.size()];//the index of Insts in each column
    //记录每一列占据的横向网格数
    vector<int> widthVec;
    int vecNum = 0;
    int totalWidth = 0;
    int tmpCol=0;
    //含有元件最多的列的元件数，即整个电路图的最大行数
    int maxLine=0;

    relocation(insts,streams,vecArr,widthVec,vecNum,totalWidth,tmpCol,maxLine);
    newWiring(vecArr,vecNum,tmpCol,totalWidth,maxLine);
    //printLines();

    outputInst(insts,"inst_out.json");
    outputStream("net_out.json");
}
