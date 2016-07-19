#include "gtest/gtest.h"
#include "../../carl/numbers/numbers.h"

TEST(GMP, Debug)
{
	carl::sint i = 4;
	mpq_class n(i);
}

TEST(GMP, Rationalize)
{
    EXPECT_TRUE( carl::rationalize<mpq_class>(carl::toDouble(mpq_class(1)/mpq_class(3))) != mpq_class(1)/mpq_class(3) );
    EXPECT_TRUE( carl::rationalize<mpq_class>(carl::toDouble(mpq_class(1)/mpq_class(20))) != mpq_class(1)/mpq_class(20) );
}
