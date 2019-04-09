class BinPos
{
private:
	friend class boost::serialization::access;
	int xbin;
	int ybin;
	int zbin;
	
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & xbin;
        ar & ybin;
        ar & zbin;
    }

public:
	BinPos() : xbin(0), ybin(0), zbin(0) {}
	BinPos(int x, int y, int z) : xbin(x), ybin(y), zbin(z){}
	~BinPos(){}

	int xBin() { return xbin; }
	int yBin() { return ybin; }
	int zBin() { return zbin; }

	void print()
	{
		std::cout << xbin << ", " << ybin << ", " << zbin << std::endl;
	}
};