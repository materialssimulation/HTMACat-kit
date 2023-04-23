from abc import ABC, abstractmethod
from ase import Atoms


class Structure(ABC):
    @abstractmethod
    def construct(self) -> Atoms:
        pass

    @abstractmethod
    def out_file_name(self) -> str:
        pass

    @abstractmethod
    def out_print(self) -> str:
        pass
