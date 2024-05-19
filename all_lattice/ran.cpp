#include <iostream>
#include <random>

int generateRandomNumber(int min, int max) {
    // 每次调用函数时都重新初始化随机数生成器
    std::mt19937 generator(std::random_device{}());

    // 创建一个均匀分布的随机数生成器
    std::uniform_int_distribution<int> distribution(min, max);

    // 生成并返回一个随机数
    return distribution(generator);
}

int main() {
    int min = 1;
    int max = 100;

    // 第一次调用
    int randomNumber1 = generateRandomNumber(min, max);
    std::cout << "First Random Number: " << randomNumber1 << std::endl;

    // 第二次调用
    int randomNumber2 = generateRandomNumber(min, max);
    std::cout << "Second Random Number: " << randomNumber2 << std::endl;

    int randomNumber3 = generateRandomNumber(min, max);
    std::cout << "3rd Random Number: " << randomNumber3<< std::endl;



    return 0;
}