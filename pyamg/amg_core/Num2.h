class Num
{
 private:
 int num;
 public:
 Num(int n);
 int getNum();
}; 

Num::Num(int n): num(n) {}
int Num::getNum()
{
 return num;
} 
