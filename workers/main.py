import asyncio


async def main():
    print("Worker stub is running. Add tasks here.")
    while True:
        await asyncio.sleep(60)


if __name__ == "__main__":
    asyncio.run(main())
