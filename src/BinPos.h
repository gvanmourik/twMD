class BinPos
{
// private:
// 	int xBin;
// 	int yBin;
// 	int zBin;

public:
	int xBin;
	int yBin;
	int zBin;

	BinPos(int x, int y, int z) : xBin(x), yBin(y), zBin(z){}
	~BinPos(){}

	void print()
	{
		std::cout << xBin << ", " << yBin << ", " << zBin << std::endl;
	}
};